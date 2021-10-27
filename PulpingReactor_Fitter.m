%% Delignification of spruce wood chips data fitting
% This file fits GPC data of DES black liquor with the pulping reactor
% model to yield kinetic data. B.J.B. Meester 2021.
clc; clear; close all;

%% Discretization GPC-chromatogram
monomer = 200;              %Lignin monomeric weight (g/mol)
DPmax = 15;                 %Maximum degree of polymerisation of lignin molecule
flowrate = 1:4;             %flow rate [ml/min]
ncount = 0;                 %loop counter
for kkkk = 1:length(flowrate)
    ncount = ncount+1;
    GPC = load([num2str(flowrate(ncount)),'T6.mat']); %Load GPC data 

    cumwd = GPC.vwd_cumdist;     %Cumulative mass distribution [%]
    mw = GPC.mw;                 %Molecular weight [g/mol]

% weight fraction of the molecular weight distribution 
% is the slope of the cum_distrubtion
    wt = abs(gradient(cumwd));  %weight fraction [%]

    st_point = 3;               %DP starting point of discretation 

    DPmaxrow = 0:DPmax+st_point;%create iteration matrix
    n = 0;
        for z = 1:DPmax
            n = n + 1;
            rowNum = find(mw>DPmaxrow(n+st_point)*monomer ... 
            & mw<DPmaxrow(n+1+st_point)*monomer); %find indices of DP slice
            DPwt(n) = sum(wt(rowNum)); %compute total weight of DP slice 
        end

    DPwt = DPwt./(sum(DPwt)); %compute weight fraction of DP slice
    DPwt(DPmax-st_point+2:DPmax) = []; %remove excess fractions
    
    fill = interp1(st_point:DPmax,DPwt,1:st_point-1,'pchip','extrap'); %interpolate low molecular weight data      
    DPwts(ncount,:) = [fill DPwt]; %add interpolated molecular weight data

end

mwp = mw; 
mwp(mwp < monomer) = NaN;   %show GPC data for M > monomer
DPslice = [1:st_point,st_point+1:DPmax];
figure(1) %plotting discretized GPC data
subplot(2,1,1)
semilogx(mwp,wt) %Plot GPC data
xlabel('Apparent molecular weight [g/mol]')
ylabel('Weight percent [%]')
title('GPC data for molecular weights higher than molecular weight monomer')
grid on
subplot(2,1,2) %Plot DP vs weight fraction
plot(DPslice,DPwts)
xlabel('Degree of polymerisation [-]')
ylabel('Weight fraction [-]')
title('Mass distribution GPC discretization')
grid on

%% Model fitting
kk = 4.1e-4;      % initial guess cracking kinetic constant       [s^-1]
kp = 1.2e-2;      % initial guess polymerisation kinetic constant [m^3 mol^-1 s^-1]
k = [kk kp]; 
D = 1e-7;         % initial guess diffusion constant        [m^2 s^-1]

km = 4.8e-13;     % mass transfer coefficient [m^2 s^-1]

DLi = 30;         % to be fitted delignification [%]

options = optimoptions('lsqnonlin','TolFun',1e-10,'StepTolerance',1e-20,'MaxIter',200,'MaxFunEvals',100); %solver options

[A,resnorm,residual,exitflag,output,lambda,jacobian] = ... 
    lsqnonlin(@(D) (funk(k,DPmax,D,km,flowrate,st_point)-DLi),D,0,1,options); %call fitter

conf_int_D = nlparci(A,residual,'jacobian',jacobian); % upper and lower 95% confidence intervals parameter D
D = A;                                                % fitted diffusion coefficient [m2 s^-1]

[A,resnorm,residual,exitflag,output,lambda,jacobian] = ... 
     lsqnonlin(@(k) (fun(k,DPmax,D,km,flowrate,st_point)-DPwts),k,[0 0],[Inf Inf],options); %call fitter kinetic rates

conf_int = nlparci(A,residual,'jacobian',jacobian); %upper and lower 95% confidence intervals parameter A(1) and A(2)
K = [A(1) A(2)]; %fitted kinetic constants

 %% Plotting solutions
[DPwt_model,t,Csolid,delignification,gin,gout,outDES,Csurface] = PulpingReactor(K(1),K(2),D,km,DPmax,flowrate,st_point); %generate mass distibution

% Determine error by inputting upper and lower kinetic rate bound
[DPwt_model_min] =  PulpingReactor(conf_int(1,1),conf_int(2,1),D,km,DPmax,flowrate,st_point); %generate lower bound mass distibution
[DPwt_model_plus] = PulpingReactor(conf_int(1,2),conf_int(2,2),D,km,DPmax,flowrate,st_point); %generate upper bound mass distibution
% Calculate error residuals
error = abs([DPwt_model(:,st_point:end)-DPwt_model_min(:,st_point:end);DPwt_model_plus(:,st_point:end)-DPwt_model(:,st_point:end)]);

ns = 0;                 % loop counter
for n = 1:length(flowrate)
    ns = ns+1;
    
    f = figure(5+ns);   % plot mass distribution
    f.Position = [593 530 706 449]; % set figure size
    plot(st_point:DPmax,DPwts(ns,st_point:end),'dk',st_point:DPmax,DPwt_model(ns,st_point:DPmax),'k')
    hold on
    shadedErrorBar(st_point:DPmax,DPwt_model(ns,st_point:end),error(ns:ns+1,:)) %plot shaded error bar
    hold off
    ylim([0 ceil((max([DPwt DPwt_model(st_point:end) DPwt_model_plus(st_point:end) DPwt_model_min(st_point:end)]))*107)/100]) %compute y-range
    ylim([0 0.15])
    formatSpec = '%.1e'; %format kinetic rate output in graph title                          
    legend('Data','Best fit')
    xlabel('Degree of polymerisation')
    ylabel('Weight fraction [-]')
    title(['Distribution DES-phase lignin with fitted parameters k_c = ' num2str(K(1),formatSpec), ' s^{-1} k_p = ' num2str(K(2),formatSpec), ' m^{3} mol ^{-1} s^{-1}']) 
    grid on
    xlim([0 DPmax])    
end

disp(['Achieved delignification = ',num2str(delignification',4),'%']) %display 

%% Fitter functions
function dif = fun(k,DPmax,D,km,phi,st_point) %function for fitting kinetic values to DP vs weight fractions
dif = PulpingReactor(k(1),k(2),D,km,DPmax,phi,st_point);
end

function DL = funk(k,DPmax,D,km,phi,st_point) %function for fitting D value to delignification percent of wood
[~,~,~,DL] = PulpingReactor(k(1),k(2),D,km,DPmax,phi,st_point);
end