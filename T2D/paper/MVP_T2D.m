function [f,sys] = MVP_T2D(par)

import casadi.*

% 4 compartment Kandarian Model modified for T2D (Aradottir et al. 2019)
% Written by: D. Krishnamoorthy, Jan 2020

% states 
I_s = MX.sym('I_s');    % Subcutaneous Insulin I_sc [U/day]
I_p = MX.sym('I_p');    % Plasma Insulin I_p [U/day]
I_e = MX.sym('I_e');    % Insulin effect on glucose I_eff [U/day]
G = MX.sym('G');        % Glucose concentration in plasma [mmol/L]

% input
u = MX.sym('u');        % exogenous insulin input [U/day]

% disturbances
SI = MX.sym('SI');      % Insulin sensitivity [1/U]
pEGP = MX.sym('pEGP');  % rate of endogenous glucose production [mmol/L day]
B = MX.sym('B');        % Endogenous insulin production co-eff beta [U L/mmol day]

% parameters
tau1 = MX.sym('tau1');  % time constant [day]
p2 = MX.sym('p2');      % delay in insulin action [1/day]
pGEZI = MX.sym('pGEZI');% rate of glucose elimination from plasma [1/day]

% ODE model
dx1 = (u - I_s)/tau1;
dx2 = (I_s - I_p)/tau1;
dx3 = p2*(I_p + B*G) - p2*I_e;
dx4 = -(pGEZI + SI*I_e)*G + pEGP;

sys.dx = vertcat(dx1,dx2,dx3,dx4);
sys.x = vertcat(I_s,I_p,I_e,G);
sys.u = u;
sys.d = vertcat(SI,pEGP,B);


sys.dx = substitute(sys.dx,tau1,par.tau1);
sys.dx = substitute(sys.dx,p2,par.p2);
sys.dx = substitute(sys.dx,pGEZI,par.pGEZI);


ode = struct('x',sys.x,'p',vertcat(sys.u,sys.d),'ode',sys.dx,'quad',0);  
opts = struct('tf',par.Ts);

f = Function('f',{sys.x,sys.u,sys.d},{sys.dx},{'x','u','d'},{'xdot'});

