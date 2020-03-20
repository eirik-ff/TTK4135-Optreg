
%% System definitions
% Discrete LTI state-space model
A = [0    0    0;
     0    0    1;
     0.1 -0.79 1.78]; 
B = [1; 0; 0.1];
C = [0 0 1];

x0 = [0;0;1];

% Number of states and inputs
mx = size(A,2);
mu = size(B,2);

% Time horizon
N = 30;

%% Objective function definitons
% Cost function matrices
Q = [0 0 0; 0 0 0; 0 0 2];
r = 1;
R = 2*r; 

% Objective function matrix
G = blkdiag(kron(eye(N),Q),kron(eye(N),R));

% Eq constraints
beq = zeros(N*mx,1);
beq(1:mx) = A*x0;

a1 = kron(eye(N),eye(mx));
a2 = kron(diag(ones(N-1,1),-1),-A);
a  = a1 + a2;
b  = kron(eye(N),-B);
Aeq = [a, b];


%% Find optimal path with quadprog
c = zeros(size(Aeq,2),1);
Aineq = zeros(size(Aeq));
bineq = zeros(size(beq));

xlb = -Inf * ones(mx,1);
xub = +Inf * ones(mx,1);
ulb =   -1 * ones(mu,1);
uub =   +1 * ones(mu,1);

vlb = [kron(ones(N,1),xlb); kron(ones(N,1),ulb)];
vub = [kron(ones(N,1),xub); kron(ones(N,1),uub)];

z0 = zeros(size(c));
z0(1:mx) = x0;

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');

[zopt,fval,exitflag,output] = quadprog(G,c,Aineq,bineq,Aeq,beq,vlb,vub,z0,options);

xopt = zopt(1:N*mx);
uopt = zopt(N*mx+1:end);

Cblk = kron(eye(N),C);
yopt = Cblk * xopt;


%% Plot quadprog solution
t = 0:N-1;


figure(1);

subplot(211);
stairs(t,uopt,'k');
grid on;
ylabel('u*');
ylimits = ylim;
ylim(ylimits + [-0.1 0.1]);
title({'TTK4135 Assignment 6 problem 3a:', ...
       'Open-loop without input blocking', ...
       sprintf('N = %i  (iterations = %i)',N,output.iterations)});

subplot(212);
plot(t,yopt);
grid on;
ylabel('y*');

xlabel('timestep');