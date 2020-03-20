%% Problem definition
% x = [x1; x2];

G = 2 * eye(2);
c = [-2; -5];

A = [1 -2;
     -1 -2;
     -1 2;
     1 0;
     0 1];
b = -[2 6 2 0 0]';

alphak = @(ai, bi, xk, pk) (bi - ai' * xk) / (ai' * pk);


%% Iteration 0
disp('----------------------------');
disp('Iteration 0');
x0 = [2 0]';
W0 = [3];
g0 = G * x0 + c;
% p0 = quadprog(G, g0, zeros(2), zeros(2,1), A(W0,:), zeros(size(W0)), -Inf * ones(2,1), +Inf * ones(2,1), x0)
plambda0 = [G, -A(W0,:)'; A(W0,:), zeros(size(A(W0,:),1))] \ [-g0; zeros(length(W0),1)];
p0 = plambda0(1:end-length(W0))
lambda0 = plambda0(end-length(W0)+1:end);

alpha0s = ones(size(b));
for i = 1:size(b,1)
    if A(i,:) * p0 < 0  % only check when ai^T pk < 0
        alpha0s(i) = alphak(A(i,:)', b(i), x0, p0);
    end
end

[alpha0, i] = min(alpha0s)


%% Iteration 1
disp('----------------------------');
disp('Iteration 1');
x1 = x0 + alpha0 * p0;
W1 = [3];
g1 = G * x1 + c;
% p1 = quadprog(G, g1, zeros(2), zeros(2,1), A(W1,:), zeros(size(W1)), -Inf * ones(2,1), +Inf * ones(2,1), x1)
% lambda_hat = linsolve(A(W1,:)', g1)

plambda1 = [G, -A(W1,:)'; A(W1,:), zeros(size(A(W1,:),1))] \ [-g1; zeros(length(W1),1)];
p1 = plambda1(1:end-length(W1))
lambda1 = plambda1(end-length(W1)+1:end)




%% Iteration 2
disp('----------------------------');
disp('Iteration 2');
x2 = x1;
W2 = [];
g2 = G * x2 + c;
% p2 = quadprog(G, g2, zeros(2), zeros(2,1), A(W2,:), zeros(size(W2)), -Inf * ones(2,1), +Inf * ones(2,1), x2)

plambda2 = [G, -A(W2,:)'; A(W2,:), zeros(size(A(W2,:),1))] \ [-g2; zeros(length(W2),1)];
p2 = plambda2(1:end-length(W2))
lambda2 = plambda2(end-length(W2)+1:end);

alpha2s = ones(size(b));
for i = 1:size(b,1)
    if A(i,:) * p2 < 0  % only check when ai^T pk < 0
        alpha2s(i) = alphak(A(i,:)', b(i), x2, p2);
    end
end

[alpha2, i] = min(alpha2s)


%% Iteration 3
disp('----------------------------');
disp('Iteration 3');
x3 = x2 + alpha2 * p2;
W3 = [1];
g3 = G * x3 + c;
% p3 = quadprog(G, g3, zeros(2), zeros(2,1), A(W3,:), zeros(size(W3)), -Inf * ones(2,1), +Inf * ones(2,1), x3)

plambda3 = [G, -A(W3,:)'; A(W3,:), zeros(size(A(W3,:),1))] \ [-g3; zeros(length(W3),1)];
p3 = plambda3(1:end-length(W3))
lambda3 = plambda3(end-length(W3)+1:end)


alpha3s = ones(size(b));
for i = 1:size(b,1)
    if A(i,:) * p3 < 0  % only check when ai^T pk < 0
        alpha3s(i) = alphak(A(i,:)', b(i), x3, p3);
    end
end

[alpha3, i] = min(alpha3s)

