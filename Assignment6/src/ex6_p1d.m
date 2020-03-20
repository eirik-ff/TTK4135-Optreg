% TTK4135 - Exercise 6, problem 1d

%% System definition
T = 0.5;
A = [1 T; 0 1];
B = [1/2*T^2; T];

Q = 2*eye(2);
R = 2;

%% Solving Riccati equation
% The Q and R matrices are halfed as dlqr's cost function doesn't include
% the 1/2 coefficient. 
[K,P,CLeigs] = dlqr(A,B,Q./2,R./2,[]);


%% Check stability
stable = all(abs(CLeigs) < 1);
if stable == 1
    disp('System is stable');
else
    disp('System is unstable');
end