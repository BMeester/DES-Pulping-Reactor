function [C_DES,wt_fraction,out] = PulpingReactor_DES(kk,kp,km,Cs,DPmax,phi,st_point,timesteps)
% This file computes the delignification of spruce wood chips in a PFR in the DES phase.
% B.J.B. Meester 2021

%% Parameters
% Constant parameters
L = 0.269;              % reactor length [m]
d = 19.86e-3;           % intrinsic reactor diameter [m]
mW = 200;               % molecular weight [g/mol]
V = L*0.25*pi*d^2;      % reactor volume [m3]
Ar = 0.25*pi*d^2;       % reactor area [m2]
% Reactor parameters
uf = phi/Ar;            % superficial velocity [m/s]

% Wood plates model
eta = 0.13;             % void fraction [-]
Lw = 0.25;              % chip length [m]
Ww = 0.03;              % chip width [m]
A = 2*Lw*Ww;            % chip area [m2]
nC = 1;                 % number of chips [-]
A = nC*A;               % total chip area [m2]
a = A/(V);              % total chip area per reactor volume [m2/m3]

% Computational parameters
n = DPmax;              % number of iterations

% Initial conditions
init(1:n) = Cs(1:n);

%% Calling solver
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'NormControl','on','Stats','off');
warning('off')
Lstep = linspace(0,L,timesteps);
[z,x] = ode45(@eq_delignin,Lstep,init,options);

% Function solution outputs
C_DES = x;                                                           % concentration DES phase lignin
wt_fraction = (x(end,:).*((1:n).*mW))./(sum(x(end,:).*((1:n).*mW))); % weight fraction lignin species in DES 
C_exit = x(end,:);                                                   % concentration outgoing lignin [mol/m3]
out = sum(C_exit.*((1:n).*mW));                                      % total outgoing lignin [g/m3]

%% Functions
    function eq = eq_delignin(y,x) 
        eq = zeros(n,1);
        nn = 0;                 
            for bb = 1:n
                nn = nn+1;
         if nn == 1
            pol_formation = 0;
            cra_consumption = 0;
        else 
            pol_formation = kp.*(sum(x(nn+1:n))).*sum(x(1:nn-1));
            cra_consumption = kk*x(nn);
        end
        if nn == n
            pol_consumption = 0;
            cra_formation = 0;
        else 
            pol_consumption = 2*kp*x(nn)*sum(x(1:n-nn));
            cra_formation = kk.*sum(x(nn+1:n));
        end
                P =  pol_formation - pol_consumption;
                Rc = cra_formation - cra_consumption;
            
                eq(nn) = (Rc + P + km*a*eta*(Cs(nn)-x(nn)))./uf;
            end      
        eq = eq(:); 
    end
end
