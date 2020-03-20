%% System definition
% Parameters
k1 = 1;
k2 = 1;
k3 = 1;

T = 0.1;

% System matrices
A = [1, T; -k2*T, 1-k1*T];
B = [0; k3*T];
C = eye(nx);

% System sizes
nx = size(A,2);
ny = size(C,1);
nu = size(B,2);

% Initial values
x0 = [5; 1];


%% Problem 2b: state observer
% State observer
Obsv_eigs = [0.5 + 0.03i, 0.5 - 0.03i];
% use A.' and C.' because this is meant to find state feedback, not observer 
% gain. See documentation 
Kf = place(A.',C.',Obsv_eigs).'; 


%% Problem 3a: Output feedback MPC
N = 50;  % timesteps
M = 10;  % time horizon in MPC

Q = diag([4, 4]);
R = 1;

u    = zeros(nu, N);
x    = zeros(nx, N+1);
y    = zeros(ny, N);

x(:,1)    = x0;

% Optimization formulation
H = blkdiag(kron(eye(M),Q),kron(eye(M),R));

% Eq constraints
beq = zeros(M*nx,1);
beq(1:nx) = A*x0;

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
    beq(1:nx) = A*x(:,t);
    
    % Calculate optimal path/input
    zopt = quadprog(H,[],[],[],Aeq,beq,vlb,vub,[],options);
    xopt = zopt(1:M*nx);
    uopt = zopt(M*nx+1:end);
    
    % Extract optimal input and apply to system
    u(:,t) = uopt(1:nu);
    x(:,t+1) = A*x(:,t) + B*u(:,t);
end  % for


%% Plot
t = 0:N;

figure(1);
set(gcf, 'position', [1000, 100, 650, 850])

subplot(311);
hold on;
grid on;
stairs(t(1:end-1),u,'LineWidth',1);
ylabel('$u$','Interpreter','latex');

title({'TTK4135 Assignment 7 problem 3b','MPC using full state feedback'});

subplot(312);
hold on;
grid on;
plot(t,x(1,:),'LineWidth',1);
ylabel('$x_1$','Interpreter','latex','FontSize',10);
legend('$x_1$','Interpreter','latex','FontSize',12);

subplot(313);
hold on;
grid on;
plot(t,x(2,:),'LineWidth',1);
ylabel('$x_2$','Interpreter','latex','FontSize',10);
legend('$x_2$','Interpreter','latex','FontSize',12,'Location','southeast');

xlabel('time step');




