%% QP definition
% x  = [A; B];

% Objective
G = [0.8 0; 0 0.4];
c = [-3; -2];

% Ineq. constraints
A = [-2 -1;
     -1 -3;
      1  0;
      0  1];
b = [-8; -15; 0; 0];

f = @(a,b) -1/2 * [a;b]' * G * [a;b] + c'*[a; b];


%% Meshgrid and contour
[X,Y] = meshgrid(0:0.1:15,0:0.1:8);
Z = zeros(size(X));
% TODO: Is there an easier way to do this?
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = f(X(i,j),Y(i,j));
    end
end

figure(1);
hold on;
title({'TTK4135 Assignment 4: Problem 2', 'Contour with constraints (in red) and iterations'});
xlabel('A [1000 kg]');
ylabel('B [1000 kg]');

contour(X,Y,Z,30);
colormap winter;
colorbar; 

xlim([0,10]);
ylim([0,6.5]);


%% Constraints
% see separate file for function: ploteqconstraints
ploteqconstraints(A, b, 'r', 'LineWidth', 1.2);

% x1 >= 0 constraint
xlimits = xlim;
lx = linspace(xlimits(1), xlimits(2), 10);
plot(lx, zeros(size(lx)), 'r', 'LineWidth', 1.2);

% x2 >= 0 constraint
ylimits = ylim;
ly = linspace(ylimits(1), ylimits(2), 10);
plot(zeros(size(ly)), ly, 'r', 'LineWidth', 1.2);



%% Plot iteration sequence
xiter = [0, 0; 
         2.4, 3.2; 
         2.25, 3.5]'
     
for i = 1:size(xiter, 2)
    p = [xiter(1,i), xiter(2,i)];
    ptxt = sprintf('%d: (%.1f,%.1f)', i, p(1), p(2));
    xshift = 0.3;
    yshift = 0.3;
    text(p(1) + xshift, p(2) + yshift, ptxt);
    
    if i == size(xiter, 2)
        plot(p(1), p(2), 'g*', 'MarkerSize', 15, 'LineWidth', 1.5);
    else
        plot(p(1), p(2), '+', 'MarkerSize', 15, 'LineWidth', 1.5);
    end
end