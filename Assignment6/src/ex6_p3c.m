
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
nublocks = [1,1,2,4,8,14];  % number of u blocks used in input blocking

%% Objective function definitons
% Cost function matrices
Q = [0 0 0; 0 0 0; 0 0 2];
r = 1;
R = 2*r; 

% Objective function matrix
G = blkdiag(kron(eye(N),Q),kron(eye(length(nublocks)),R));

% Eq constraints
beq = zeros(N*nx,1);
beq(1:nx) = A*x0;

a1 = kron(eye(N),eye(nx));  % diagonal identity matrices
a2 = kron(diag(ones(N-1,1),-1),-A);  % off-diagonal A matrices
a  = a1 + a2;

b = zeros(N,length(nublocks));
t = 0;
for i = 1:length(nublocks)
    il = t + 1;
    t = t + nublocks(i);
    iu = t;
    
    b(il:iu,i) = ones(nublocks(i),1);
end
b  = kron(b,-B);

Aeq = [a, b];



%% Find optimal path with quadprog
xlb = -Inf * ones(nx,1);
xub = +Inf * ones(nx,1);
ulb =   -1 * ones(nu,1);
uub =   +1 * ones(nu,1);

vlb = [kron(ones(N,1),xlb); kron(ones(length(nublocks),1),ulb)];
vub = [kron(ones(N,1),xub); kron(ones(length(nublocks),1),uub)];

z0 = zeros(size(Aeq,2),1);
z0(1:nx) = x0;

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');

[zopt,fval,exitflag,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,z0,options);

xopt = zopt(1:N*nx);
u    = zopt(N*nx+1:end);
uopt = zeros(N,1);
t = 0;
for i = 1:length(nublocks)
    il = t + 1;
    t = t + nublocks(i);
    iu = t;
    
    uopt(il:iu) = u(i) * ones(nublocks(i),1);
end

Cblk = kron(eye(N),C);
yopt = Cblk * xopt;



%% Plot quadprog solution
t = 0:N-1;


figure(1);

subplot(211);
stairs(t,uopt,'k');
grid on;
hold on;
for i = nublocks
    xline(i,'--');
end
ylabel('u*');
ylimits = ylim;
ylim(ylimits + [-0.1 0.1]);
title({'TTK4135 Assignment 6 problem 3c:', ... 
       'Input blocking with increasing step length (shown by dashed lines)', ...
       sprintf('N = %i  (iterations = %i)',N,output.iterations)});

subplot(212);
plot(t,yopt);
grid on;
ylabel('y*');

xlabel('timestep');