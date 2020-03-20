%% LP definition
% x  = [A; B; x3; x4];
% x3,x4: slack variables
Aeq = [2 1 1 0; 
       1 3 0 1];
beq = [8; 15];
c = [-3/2; -1; 0; 0];

f = @(a,b) c'*[a; b; 0; 0];


%% Meshgrid and contour
[X,Y] = meshgrid(0:1:15,0:1:8);
Z = zeros(size(X));
% TODO: Is there an easier way to do this?
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = f(X(i,j),Y(i,j));
    end
end

figure(1);
hold on;
title('TTK4135 Assignment 3: Problem 2');
xlabel('A [1000 kg]');
ylabel('B [1000 kg]');

contour(X,Y,Z,15);
colorbar; 


%% Constraints
% see separate file for function: ploteqconstraints
ploteqconstraints(Aeq, beq, 'r', 'LineWidth', 1.2);


%% Simplex method
x0 = [0 0 beq']';
[xopt, fval, iterations] = simplex(c, Aeq, beq, x0, 'report');


%% Plot x* and iterations
for i = 1:size(iterations, 2)
    p = [iterations(1,i), iterations(2,i)];
    ptxt = sprintf('%d: (%.1f,%.1f)', i, p(1), p(2));
    xshift = 0.3;
    yshift = 0.3;
    text(p(1) + xshift, p(2) + yshift, ptxt);
    
    if i == size(iterations, 2)
        plot(p(1), p(2), 'g*', 'MarkerSize', 15, 'LineWidth', 1.5);
    else
        plot(p(1), p(2), '+', 'MarkerSize', 15, 'LineWidth', 1.5);
    end
end

