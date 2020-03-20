%% Rosenbrock
x = sym('x', [2,1]);
f_sym = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
G_sym = gradient(f_sym, x);
H_sym = hessian(f_sym, x);

f = matlabFunction(f_sym, 'vars', {x});
G = matlabFunction(G_sym, 'vars', {x});
H = matlabFunction(H_sym, 'vars', {x});

x0_1 = [1.2; 1.2];
x0_2 = [-1.2; 1];

%% Compute steepest descent
x0 = x0_2;
tic;
[x_opt, fval_opt, x_iter, f_iter, alpha] = SteepestDescent(x0,f,G);
disp('Steepest descent:');
toc
tic;
[x_opt_n, fval_opt_n, x_iter_n, f_iter_n, alpha_n] = NewtonsDescent(x0,f,G,H);
disp('Newton:');
toc


%% Plot steepest
figure(1);
clf;
hold on;
subplot(311);
plot_iter_rosenbrock(x_iter);
title({'Optreg Assignment 9 problem 1a', ...
       'Iterations of steepest descent on the Rosenbrock function',...
       sprintf('using initial point [%.2f %.2f]^T', x0(1), x0(2))});

subplot(312);
semilogy(1:length(alpha), alpha);
hold on;
grid on;
title('Step length \alpha_k');
xlabel('iteration k');
ylabel('step length \alpha');


subplot(313);
semilogy(1:length(f_iter), f_iter);
hold on;
grid on;
title('Function value f(x_k)');
xlabel('iteration k');
ylabel('f(x_k)');


%% Plot newtons
figure(2);
clf;
hold on;
subplot(311);
plot_iter_rosenbrock(x_iter_n);
title({'Optreg Assignment 9 problem 1a', ...
       'Iterations of Newtons descent on the Rosenbrock function',...
       sprintf('using initial point [%.2f %.2f]^T', x0(1), x0(2))});

subplot(312);
semilogy(1:length(alpha_n), alpha_n);
hold on;
grid on;
title('Step length \alpha_k');
xlabel('iteration k');
ylabel('step length \alpha');


subplot(313);
semilogy(1:length(f_iter_n), f_iter_n);
hold on;
grid on;
title('Function value f(x_k)');
xlabel('iteration k');
ylabel('f(x_k)');
