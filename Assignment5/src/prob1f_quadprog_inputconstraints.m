% FILE: prob1f_quadprog_inputconstraints.m

run prob1_defs


%% Find optimal path with quadprog
options = optimoptions(@quadprog,'Algorithm','interior-point-convex');

[zopt,fval,exitflag,output] = quadprog(G,c,Aineq,bineq,Aeq,beq,vlb,vub,z0,options);

xopt = zopt(1:N*mx);
uopt = zopt(N*mx+1:end);

Cblk = kron(eye(N),C);
yopt = Cblk * xopt;


%% Plot quadprog solution
t = 0:N-1;


figure(1);

subplot(211);
stairs(t,uopt,'k');
grid on;
ylabel('u*');
title({'TTK4135 Assignment 5 problem 1f: Solution using MATLABs quadprog with input constraints', 
       sprintf('r = %G  (algorithm: "%s", iterations = %i)',r,output.algorithm,output.iterations)});

subplot(212);
plot(t,yopt);
grid on;
ylabel('y*');

xlabel('timestep');