% FILE: prob1d_kkt_sol.m

run prob1_defs


%% Solve KKT system (16.4)
H = [G, -Aeq'; Aeq, zeros(min(size(Aeq)))];
h = [zeros(size(Aeq,2),1); beq];

zlambda = H \ h;
zopt      = zlambda(1:N*mx+N*mu);
lambdaopt = zlambda(N*mx+N*mu:end);

xopt = zopt(1:N*mx);
uopt = zopt(N*mx+1:end);

Cblk = kron(eye(N),C);
yopt = Cblk * xopt;

%% Plot KKT system solution
t = 0:N-1;

figure(1);

subplot(211);
stairs(t,uopt,'k');
grid on;
ylabel('u*');
title({'TTK4135 Assignment 5 problem 1d: Solution to KKT system (16.4)','Optimal control sequence and optimal trajectory in open-loop'});

subplot(212);
plot(t,yopt);
grid on;
ylabel('y*');

xlabel('timestep');