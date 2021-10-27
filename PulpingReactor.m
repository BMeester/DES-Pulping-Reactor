function [wtf,tt,Csolid,delignification,gin,gout,outDES,Csurface] = PulpingReactor(kk,kp,D,km,DPmax,phi,st_point)
% This file computes the delignification of spruce wood chips in a PFR.
% B.J.B. Meester 2021

%% Parameters
% Constant parameters
L = 0.269;          % reactor length [m]
d = 19.86e-3;       % intrinsic reactor diameter [m]
mW = 200;           % monomer molecular weight [g/mol]
V = L*0.25*pi*d^2;  % reactor volume [m3]
Ar = 0.25*pi*d^2;   % reactor area [m2]
t = 180;            % reaction time [min]
t_tot = t*60;       % reaction time [s]
Woodin = 10;        % spruce wood in reactor at time = 0 [g]

% Reactor parameters
phi = phi./60*10^-6;    % flowrate [m3/s] (input= ml/min*conversion factor)
uf = phi./Ar;           % superficial velocity [m/s]
t_res = V./phi;         % residence time [s]

% Wood Plates model
eta = 0.13;             % void fraction [-]
Lw = 0.25;              % chip length [m]
Ww = 0.03;              % chip width [m]
A = 2*Lw*Ww;            % chip area [m2]
nC = 1;                 % number of chips [-]
A = nC*A;               % total chip area [m2]
a = A/(V);              % total chip area per reactor volume [m2/m3]

% Computational parameters
n = DPmax;              %number of iterations
Cwood = [7.0723E-05,0.026370495,0.058865064,0.035362575,0.030122451,0.039848341,0.049896177,0.046265472,0.038854025,0.033614917,0.029743237,0.028446201,0.02835567,0.025450231,0.025198676,0.024787421,0.022765525,0.023319698,0.020896312,0.02145079,0.020555345,0.020706502,0.017799944,0.019721544,0.018345165,0.017128183,0.017440261,0.015580711,0.015409007,0.015641031,0.015550805,0.014092847,0.014002723,0.012303688,0.013903139,0.012123338,0.011952346,0.010333466,0.011772301,0.010153411,0.010063582,0.009732014,0.009561724,0.009390722,0.008014293,0.008969325,0.007592174,0.007179799,0.00708997,0.008206636]./1000;  % Initial concentration of lignin (Ascending DP) [mol/L]
timesteps = 100;        % solver timesteps
% Initial conditions
init = zeros(2*n,1);
init(st_point:n) = Cwood(st_point:n).*Woodin.*0.5; % internal concentration internal lignin
init(n+st_point:2*n) = 0;                          % initial concentration surface lignin

%% Calling solver
warning('off')
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'NormControl','on','Stats','off'); %solver 
wtf = zeros(length(phi),DPmax);         % weight fraction pre-allocating 
delignification = zeros(length(phi),1); % delignification pre-allocating 
nc = 0;                                 % index counter
for ind = 1:length(phi)
    nc = nc + 1;
    [z,x] = ode45(@eq_delignin,linspace(0,t_tot,timesteps),init,options);

    
    gin = sum(init(1:n).*((1:n)'.*mW));     % lignin in wood [g] at t_0
    gout = sum(x(end,1:DPmax).*((1:n).*mW));% lignin in wood [g] at t_end

    % Function output
    [~,wtf(nc,:),outDES] = PulpingReactor_DES(kk,kp,km,x(end,DPmax+1:end),DPmax,phi(nc),st_point,timesteps);
    tt = z; %time profile [s]
    Csolid = x(:,1:DPmax).*(((1:n)'.*mW)');         % concentration profile internal lignin [g/m3]
    Csurface = x(:,DPmax+1:end).*(((1:n)'.*mW)');   % concentration profile surface lignin  [g/m3]

    disp('Solution found')                          % display "solution found" to indicate calculation progess
    delignification(nc) = (gin - gout)/gin *100;    % delignification

end

%% Functions
    function eq = eq_delignin(y,x)
               % Notation:
               % x(1:DPmax) =     mole DP lignin wood internal / m3
               % x(DPmax+1:end) = mole DP lignin wood surface / m3
               eq = zeros(2*DPmax,1); % pre-allocation 
               %% Calculate DES lignin concentration
               xCDES = (PulpingReactor_DES(kk,kp,km,x(DPmax+1:end),DPmax,phi(nc),st_point,timesteps))';
               xCDES = xCDES(:,2); %take time equivalent wrt residence time  

               %% Flux equations for internal and surface lignin
               eq(1:DPmax) = -(D*Ww*Lw/V).*x(1:DPmax);  
               eq(DPmax+1:end) = (D*Ww*Lw/V).*x(1:DPmax)-(km*a*eta*(x(DPmax+1:end)-xCDES(1:DPmax))); 
                              
                eq = eq(:);
    end            
end
