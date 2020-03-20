%% System definition
% Parameters
k1 = 1;
k2 = 1;
k3 = 1;

T = 0.1;

% System matrices
A = [1, T; -k2*T, 1-k1*T];
nx = size(A,2);

B = [0; k3*T];
nu = size(B,2);

C = eye(nx);
ny = size(C,1);


% Initial values
x0 = [5; 1];


%% Problem 2b: state observer
% State observer
Obsv_eigs = [0.5 + 0.03i, 0.5 - 0.03i];
% use A.' and C.' because this is meant to find state feedback, not observer 
% gain. See documentation 
Kf = place(A.',C.',Obsv_eigs).'; 


%% Problem 4a: Riccati matrix
Q = diag([4, 4]);
R = 1;
[K,P] = dlqr(A,B,Q,R,[]);


%% Problem 3a: Output feedback MPC
N = 50;  % timesteps
M = 10;  % time horizon in MPC

u_fin = zeros(nu, N);
x_fin = zeros(nx, N+1);

u_inf = zeros(nu, N);
x_inf = zeros(nx, N+1);

x_fin(:,1) = x0;
x_inf(:,1) = x0;

% Optimization formulation
Gfin = blkdiag(kron(eye(M),Q),kron(eye(M),R)); % finite horizon
Ginf = Gfin; % infinite horizon
Ginf(M*nx-1:M*nx, M*nx-1:M*nx) = P;

% Eq constraints
beq_fin = zeros(M*nx,1);
beq_fin(1:nx) = A*x0;

beq_inf = zeros(M*nx,1);
beq_inf(1:nx) = A*x0;

a1 = kron(eye(M),eye(nx));
a2 = kron(diag(ones(M-1,1),-1),-A);
a  = a1 + a2;
b  = kron(eye(M),-B);
Aeq = [a, b];


% Lower and upper bounds
xlb = -Inf * ones(nx,1);
xub = +Inf * ones(nx,1);
ulb =   -4 * ones(nu,1);
uub =   +4 * ones(nu,1);

vlb = [kron(ones(M,1),xlb); kron(ones(M,1),ulb)];
vub = [kron(ones(M,1),xub); kron(ones(M,1),uub)];


% MPC simulation
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');
for t = 1:N
    % Use previous timestep as initial value for MPC
    beq_fin(1:nx) = A*x_fin(:,t);
    beq_inf(1:nx) = A*x_inf(:,t);
    
    % Calculate optimal path/input
    zopt_fin = quadprog(Gfin,[],[],[],Aeq,beq_fin,vlb,vub,[],options);
    xopt_fin = zopt_fin(1:M*nx);
    uopt_fin = zopt_fin(M*nx+1:end);
    
    zopt_inf = quadprog(Ginf,[],[],[],Aeq,beq_inf,vlb,vub,[],options);
    xopt_inf = zopt_inf(1:M*nx);
    uopt_inf = zopt_inf(M*nx+1:end);
    
    % Extract optimal input and apply to system
    u_fin(:,t) = uopt_fin(1:nu);
    x_fin(:,t+1) = A*x_fin(:,t) + B*u_fin(:,t);
    
    u_inf(:,t) = uopt_inf(1:nu);
    x_inf(:,t+1) = A*x_inf(:,t) + B*u_inf(:,t);
end  % for


%% Plot compare finite vs infinite horizon
t = 0:N;

figure(1);
% set(gcf, 'position', [1000, 100, 650, 850])

subplot(311);
hold on;
grid on;
stairs(t(1:end-1),u_fin,'LineWidth',1);
stairs(t(1:end-1),u_inf,'LineWidth',1);
ylabel('$u$','Interpreter','latex');

title({'TTK4135 Assignment 7 problem 4b',...
       'Infinite horizon MPC with full state feedback',...
       sprintf('vs finite horizon (problem 3b). N = %G', M)});

subplot(312);
hold on;
grid on;
plot(t,x_fin(1,:),'LineWidth',1);
plot(t,x_inf(1,:),'LineWidth',1);
ylabel('$x_1$','Interpreter','latex','FontSize',10);
legend('Finite horizon: $x_1$','Infinite horizon: $x_1$',...
       'Interpreter','latex','FontSize',12);

subplot(313);
hold on;
grid on;
plot(t,x_fin(2,:),'LineWidth',1);
plot(t,x_inf(2,:),'LineWidth',1);
ylabel('$x_2$','Interpreter','latex','FontSize',10);
legend('Finite horizon: $x_2$','Infinite horizon: $x_2$',...
       'Interpreter','latex','FontSize',12,'Location','southeast');

xlabel('time step');


