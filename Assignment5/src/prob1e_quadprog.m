% FILE: prob1e_quadprog.m

run prob1_defs

legend_list = {};
l_idx = 1;
for r = [0.1, 1, 10]
    %% Find optimal path with quadprog
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex','Display','off');
    
    % Cost function matrices
    Q = [0 0 0; 0 0 0; 0 0 2];
    R = 2*r; 

    % Objective function matrix
    G = blkdiag(kron(eye(N),Q),kron(eye(N),R));

    [zopt,~,~,output] = quadprog(G,[],[],[],Aeq,beq,[],[],z0,options);

    xopt = zopt(1:N*mx);
    uopt = zopt(N*mx+1:end);

    Cblk = kron(eye(N),C);
    yopt = Cblk * xopt;


    %% Plot quadprog solution
    t = 0:N-1;


    figure(1);

    subplot(211);
    hold on;
    stairs(t,uopt);
    grid on;
    ylabel('u*');

    subplot(212);
    hold on;
    plot(t,yopt);
    grid on;
    ylabel('y*');

    xlabel('timestep');
    
    
    legend_list{l_idx} = sprintf('r = %G (iterations = %i)',r,output.iterations);
    l_idx = l_idx + 1;
end  % for
subplot(211);
title({'TTK4135 Assignment 5 problem 1e: Solution using MATLABs quadprog', ...
       'Optimal control sequence and optimal trajectory in open-loop', ...
       sprintf('(algorithm: "%s")',output.algorithm)});
subplot(212);
legend(legend_list);