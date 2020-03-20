%% Initial points
x0_1 = [1.2; 1.2];
x0_2 = [-1.2; 1];


%% Run and plot for x0_1
figure(1);
clf;

subplot(211);
tic
[x1,fval1,x1_iter] = nelder_mead(x0_1,[]);%,'report');
disp('Nelder-Mead x0_1:');
toc
title(sprintf('Nelder-Mead on the Rosenbrock function with x_0 = [%.1f, %.1f]^T', x0_1(1), x0_1(2)));

subplot(212);
plot_iter_rosenbrock(x1_iter);


%% Run and plot for x0_2
figure(2);
clf;

subplot(211);
tic
[x2,fval2,x2_iter] = nelder_mead(x0_2,[]);%,'report');
disp('Nelder-Mead x0_2:');
toc
title(sprintf('Nelder-Mead on the Rosenbrock function with x_0 = [%.1f, %.1f]^T', x0_2(1), x0_2(2)));

subplot(212);
plot_iter_rosenbrock(x2_iter);
