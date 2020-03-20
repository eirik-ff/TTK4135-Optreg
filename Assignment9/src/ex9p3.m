%% Finite difference gradient
Gf = @(x, e) [-200*(x(2) - x(1)) - 2*(1-x(1)) + 101*e; 
              200*(x(2) - x(1)) + 100*e];
          
%% Analytic gradient
Df = @(x) [-200*(x(2) - x(1)) - 2*(1-x(1)); 
           200*(x(2) - x(1))];
       
       
%% b: calculations
x1 = [0.5; 0.5];
x2 = [1; 1];
eps = [1e-16; 1e-12; 1e-8; 1e-4; 1e-0];
% eps = logspace(-12, -1,12);

Df1  = Df(x1);
Df2  = Df(x2);

Gf1s = NaN(2,length(eps));
Gf2s = NaN(2,length(eps));

for i = 1:length(eps)
    e = eps(i);
    Gf1s(:,i) = Gf(x1, e);
    Gf2s(:,i) = Gf(x2, e);
end

% difference between finite difference and analytic
d1 = Gf1s - Df1;
d2 = Gf2s - Df2;


%% Plot difference
figure(1);
clf;

semilogy(1:length(eps), d1(1,:));
xticks('manual');
xticks(1:length(eps));
xticklabels(eps);
hold on; 
grid on;

title({'Optreg Assignment 9 problem 3', ...
       'Difference between finite difference gradient',...
       'and analytic gradient'});
xlabel('\epsilon');
ylabel('difference');

   
   


