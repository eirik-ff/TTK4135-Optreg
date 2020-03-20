
%% System definitions
% Discrete LTI state-space model
A = [0    0    0;
     0    0    1;
     0.1 -0.79 1.78]; 
B = [1; 0; 0.1];
C = [0 0 1];

x0 = [0;0;1];

% Number of states and inputs
nx = size(A,2);
nu = size(B,2);

% Time horizon
N = 30;

% Input blocking 
nublocks = 6;  % number of u blocks used in input blocking

%% Objective function definitons
% Cost function matrices
Q = [0 0 0; 0 0 0; 0 0 2];
r = 1;
R = 2*r; 

% Objective function matrix
G = blkdiag(kron(eye(N),Q),kron(eye(nublocks),R));

% Eq constraints
beq = zeros(N*nx,1);
beq(1:nx) = A*x0;

a1 = kron(eye(N),eye(nx));  % diagonal identity matrices
a2 = kron(diag(ones(N-1,1),-1),-A);  % off-diagonal A matrices
a  = a1 + a2;

b  = kron(eye(nublocks),ones(N/nublocks,1));
b  = kron(b,-B);

Aeq = [a, b];



%% Find optimal path with quadprog
xlb = -Inf * ones(nx,1);
xub = +Inf * ones(nx,1);
ulb =   -1 * ones(nu,1);
uub =   +1 * ones(nu,1);

vlb = [kron(ones(N,1),xlb); kron(ones(nublocks,1),ulb)];
vub = [kron(ones(N,1),xub); kron(ones(nublocks,1),uub)];

z0 = zeros(size(c));
z0(1:nx) = x0;

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');

[zopt,fval,exitflag,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,z0,options);

xopt = zopt(1:N*nx);
uopt = zopt(N*nx+1:end);
uopt = repelem(uopt,N/nublocks);

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
title({'TTK4135 Assignment 6 problem 3b:', ... 
       'Input blocking with 6 equisized blocks', ...
       sprintf('N = %i  (iterations = %i)',N,output.iterations)});

subplot(212);
plot(t,yopt);
grid on;
ylabel('y*');

xlabel('timestep');