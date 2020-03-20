% FILE: prob2c_mpc_two_models.m

run prob1_defs

%% Actual plan model definition
A2 = [0 0 0; 0 0 1; 0.1 -0.855 1.85];
B2 = [1; 0; 0];
C2 = [0 0 1];

%% MPC simulation
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');
Cblk = kron(eye(N),C);

xtraj = zeros(mx,N);
utraj = zeros(mu,N);
ytraj = zeros(1, N);

xtraj2 = zeros(mx,N);
utraj2 = zeros(mu,N);
ytraj2 = zeros(1, N);


zopt = z0;
x = x0;
x2 = x0;
for i = 1:N
    % IMPORTANT! Update eq constraints for each iteration to include inital
    % conditions 
    beq(1:mx) = A*x;
    
    [zopt,~,~,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,[],options);
    
    xopt = zopt(1:N*mx);
    uopt = zopt(N*mx+1:end);
    yopt = Cblk * xopt;
    
    % Simulated model
    xtraj(:,i) = A*x + B*uopt(1:mu);
    utraj(:,i) = uopt(1:mu);
    ytraj(:,i) = C*xtraj(:,i);
    
    x = xtraj(:,i); 
    
    
    
    % IMPORTANT! Update eq constraints for each iteration to include inital
    % conditions 
    beq(1:mx) = A2*x2;
    
    [zopt2,~,~,output] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,[],options);
    
    xopt2 = zopt2(1:N*mx);
    uopt2 = zopt2(N*mx+1:end);
    yopt2 = Cblk * xopt2;
    
    % Actual model
    % x_i = A*x_i-1 + B*u_i-1
    xtraj2(:,i) = A2*x2 + B2*uopt2(1:mu);
    utraj2(:,i) = uopt2(1:mu);
    ytraj2(:,i) = C2*xtraj2(:,i);
    
    x2 = xtraj2(:,i);  
end  % for


%% Plot MPC solution
t = 0:N-1;


figure(1);

subplot(211);
grid on;
hold on;
stairs(t,utraj,'r--');
stairs(t,utraj2,'b');
ylabel('u*');
title({'TTK4135 Assignment 5 problem 2c: MPC with two different models', ...
       sprintf('r = %i  (algorithm: "%s", iterations = %i)',r,output.algorithm,output.iterations)});

subplot(212);
grid on;
hold on;
plot(t,ytraj, 'r--');
plot(t,ytraj2, 'b');
ylabel('y*');
legend('Model (1)', 'Model (7)');

xlabel('timestep');

