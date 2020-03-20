%% System definition
% Parameters
k1 = 1;
k2 = 1;
k3 = 1;

T = 0.1;

% System matrices
A = [1, T; -k2*T, 1-k1*T];
B = [0; k3*T];
C = [1, 0];

% System sizes
nx = size(A,2);
ny = size(C,1);
nu = size(B,2);

% Initial values
x0 = [5; 1];
x0_hat = [6; 0];


%% Problem 2a: feedback matrix
Q = diag([4, 4]);
R = 1;
[K,P,CLeigs] = dlqr(A,B,Q/2,R/2,[]);

Acl = A - B*K;


%% Problem 2b: state observer
% State observer
Obsv_eigs = [0.5 + 0.03i, 0.5 - 0.03i];
% use A.' and C.' because this is meant to find state feedback, not observer 
% gain. See documentation 
Kf = place(A.',C.',Obsv_eigs).'; 

Aobsv = A - Kf*C;

% Full matrix phi
phi = [Acl, B*K; zeros(nx), Aobsv];
phi_eigs = eigs(phi);

% Simulation
N = 50;

u    = zeros(nu, N);
x    = zeros(nx, N+1);
y    = zeros(ny, N);
xhat = zeros(nx, N+1);

x(:,1)    = x0;
xhat(:,1) = x0_hat;

for t = 1:N
    u(:,t)   = -K*xhat(:,t);
    x(:,t+1) = A*x(:,t) + B*u(:,t);
    y(:,t)   = C*x(:,t);
    
    xhat(:,t+1) = A*xhat(:,t) + B*u(:,t) + Kf*(y(:,t) - C*xhat(:,t));
end


%% Plot
t = 0:N;

figure(1);
set(gcf, 'position', [1000, 100, 650, 850])

subplot(311);
hold on;
grid on;
stairs(t(1:end-1),u,'LineWidth',1);
ylabel('$u$','Interpreter','latex');

title({'TTK4135 Assignment 7 problem 2b','Feedback using LQR and state-observer'});

subplot(312);
hold on;
grid on;
plot(t,x(1,:),'LineWidth',1);
plot(t,xhat(1,:),'--','LineWidth',1);
ylabel('$x_1$','Interpreter','latex','FontSize',10);
legend('$x_1$','$\hat{x}_1$','Interpreter','latex','FontSize',12);

subplot(313);
hold on;
grid on;
plot(t,x(2,:),'LineWidth',1);
plot(t,xhat(2,:),'--','LineWidth',1);
ylabel('$x_2$','Interpreter','latex','FontSize',10);
legend('$x_2$','$\hat{x}_2$','Interpreter','latex','FontSize',12,'Location','southeast');

xlabel('time step');




