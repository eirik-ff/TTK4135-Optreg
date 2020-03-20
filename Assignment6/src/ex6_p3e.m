
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

a1 = kron(eye(N),eye(nx));
a2 = kron(diag(ones(N-1,1),-1),-A);
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


% Lower and upper bounds
xlb = -Inf * ones(nx,1);
xub = +Inf * ones(nx,1);
ulb =   -1 * ones(nu,1);
uub =   +1 * ones(nu,1);

vlb = [kron(ones(N,1),xlb); kron(ones(length(nublocks),1),ulb)];
vub = [kron(ones(N,1),xub); kron(ones(length(nublocks),1),uub)];

% Initial state
z0 = zeros(size(Aeq,2),1);
z0(1:nx) = x0;


options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');
Cblk = kron(eye(N),C);

%% MPC simulation
xtraj = zeros(nx,N);
utraj = zeros(nu,N);
ytraj = zeros(1, N);

uopt_openloop = [];
yopt_openloop = [];

x = x0;
for i = 1:N
    % IMPORTANT! Update eq constraints for each iteration to include inital
    % conditions 
    beq(1:nx) = A*x;
    
    [zopt,~,~,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,[],options);
    
    xopt = zopt(1:N*nx);
    uopt = zopt(N*nx+1:end);
    yopt = Cblk * xopt;
    if i == 1
        u = zeros(N,1);
        t = 0;
        for j = 1:length(nublocks)
            il = t + 1;
            t = t + nublocks(j);
            iu = t;

            u(il:iu) = uopt(j) * ones(nublocks(j),1);
        end
        
        uopt_openloop = u;
        yopt_openloop = yopt;
    end  % if
    
    % x_i = A*x_i-1 + B*u_i-1
    xtraj(:,i) = A*x + B*uopt(1:nu);
    utraj(:,i) = uopt(1:nu);
    ytraj(:,i) = C*xtraj(:,i);
    
    x = xtraj(:,i);
end  % for


%% Plot MPC solution
t = 0:N-1;


figure(1);

subplot(211);
grid on;
hold on;
stairs(t,utraj,'b');
stairs(t,uopt_openloop, 'r--');
ylabel('u*');
title({'TTK4135 Assignment 6 problem 3e:', ...
       'MPC with input blocking of increasing size', ...
       sprintf('N = %i',N)});

subplot(212);
grid on;
hold on;
plot(t,ytraj, 'b');
plot(t,yopt_openloop, 'r--');
ylabel('y*');
legend('MPC', 'Open-loop');

xlabel('timestep');

