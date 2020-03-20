
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


% Lower and upper bounds
xlb = -Inf * ones(mx,1);
xub = +Inf * ones(mx,1);
ulb =   -1 * ones(mu,1);
uub =   +1 * ones(mu,1);

vlb = [kron(ones(N,1),xlb); kron(ones(N,1),ulb)];
vub = [kron(ones(N,1),xub); kron(ones(N,1),uub)];

% Initial state
z0 = zeros(size(Aeq,2),1);
z0(1:mx) = x0;


options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');
Cblk = kron(eye(N),C);

%% MPC simulation
xtraj = zeros(mx,N);
utraj = zeros(mu,N);
ytraj = zeros(1, N);

uopt_openloop = [];
yopt_openloop = [];

x = x0;
for i = 1:N
    % IMPORTANT! Update eq constraints for each iteration to include inital
    % conditions 
    beq(1:mx) = A*x;
    
    [zopt,~,~,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,[],options);
    
    xopt = zopt(1:N*mx);
    uopt = zopt(N*mx+1:end);
    yopt = Cblk * xopt;
    if i == 1
        uopt_openloop = uopt;
        yopt_openloop = yopt;
    end  % if
    
    % x_i = A*x_i-1 + B*u_i-1
    xtraj(:,i) = A*x + B*uopt(1:mu);
    utraj(:,i) = uopt(1:mu);
    ytraj(:,i) = C*xtraj(:,i);
    
    x = xtraj(:,i);
end  % for


%% Plot MPC solution
t = 0:N-1;


figure(1);

subplot(211);
grid on;
hold on;
stairs(t,utraj,'k');
stairs(t,uopt_openloop, 'r--');
ylabel('u*');
title({'TTK4135 Assignment 6 problem 3d:', ...
       'MPC with no input blocking', ...
       sprintf('N = %i',N)});

subplot(212);
grid on;
hold on;
plot(t,ytraj, 'b');
plot(t,yopt_openloop, 'r--');
ylabel('y*');
legend('MPC', 'Open-loop');

xlabel('timestep');

