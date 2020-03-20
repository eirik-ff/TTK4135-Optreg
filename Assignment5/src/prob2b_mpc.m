% FILE: prob2b_mpc.m

run prob1_defs

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
title({'TTK4135 Assignment 5 problem 2b: MPC', ...
       sprintf('r = %i  (algorithm: "%s", iterations = %i)',r,output.algorithm,output.iterations)});

subplot(212);
grid on;
hold on;
plot(t,ytraj, 'b');
plot(t,yopt_openloop, 'r--');
ylabel('y*');
legend('MPC', 'Open-loop');

xlabel('timestep');

