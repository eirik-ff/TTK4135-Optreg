%% Definitions
Aeq = [3 2 1; 2 2 3];
beq = [7200; 6000];
c   = [-100; -75; -55];
lb  = zeros(3, 1);

% For testing sensitivity
beq_sensA = beq + [1; 0];
beq_sensB = beq + [0; 1];


%% Solve with linprog
[x, fval, exitflag, output, lambda] = linprog(c, [], [], Aeq, beq, lb);
[x_sensA, fval_sensA] = linprog(c, [], [], Aeq, beq_sensA, lb);
[x_sensB, fval_sensB] = linprog(c, [], [], Aeq, beq_sensB, lb);

display('No shift:');
x, fval

display('---------------');
display('Shift to A:');
x_sensA, fval_sensA

display('---------------');
display('Shift to B:');
x_sensB, fval_sensB

