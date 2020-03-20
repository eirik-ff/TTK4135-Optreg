function [x_opt, fval_opt, x_iter, f_iter, alpha] = steepest_descent_suggestion(x0)
% Returns 
% x_opt: x^*
% fval_opt: f(x^*)
% x_iter: all iterates k -- column k contains x_k
% f_iter: a vector of all function values f(x_k)
% alpha: a vector of all step lenghts alpha_k

% Termination criteria
maxiter = 'largest allowed number of iterations k -- some large number';
grad_tol = 'how small we want ||nabla f|| to be (close to the solution) -- some small number';

x0 = x0(:);
n = size(x0,1); % Number of variables

% Declare some variables
x =     NaN(n,maxiter);
p =     NaN(n,maxiter);
grad =  NaN(n,maxiter);
alpha = NaN(1,maxiter);
fval =  NaN(1,maxiter);

%
% Do some calcualtions before the while loop. I.e., do the first iteration:
k = 1; % iteration number
x(:,k) = 'store the initial point.';
%look at the different functions below and finish them before continuing.
fval(k) = 'store the corresponding function value';
grad(:,k) = 'store the corresponding function value'; 
p(:,k) = 'find the steepest descent direction.'; 

alpha_0 = 'such that ||alpha0 pk|| = 1';

alpha(k) = 'do the line search.';
x(:,k+1) = 'take the step';
grad(:,k+1) = 'store the gradient';
k = 'increase the counter';

while 'maxiter not exceeded' && 'the norm of the gradient larger than grad_tol'
    fval(k) = f(x(:,k));    % Evaluate the Rosenbrock function
    p(:,k) = sd(grad(:,k)); % Calculate steepest-descent direction based on gradient
    alpha_0 = 'initial guess of step lenght, see p. 59 in N&W'; 
    alpha(k) = linesearch(x(:,k), p(:,k), fval(k), grad(:,k), alpha_0); % Determine alpha using Alg. 3.1
    x(:,k+1) = '';
    grad(:,k+1) = gradient(''); 
    k = k+1;
end
fval(k) = f(x(:,k)); % Final function value

% Delete unused space
x = x(:,1:k);
p = p(:,1:k);
grad = grad(:,1:k);
alpha = alpha(1:k);
fval = fval(1:k);

% Return values
x_opt = x(:,end);
fval_opt = f(x_opt);
x_iter = x;
f_iter = fval;

end

% Function returning the steepest-descent direction based on the gradient
% of f
function p = sd(grad)
    p = 'excpression for p';
end

% Function implementing Algorithm 3.1, page 37 in N&W
function alpha_k = linesearch(xk, pk, fk, gradk, alpha_0)
    alpha = alpha_0;
    rho = 'contraction factor';
    c1 = 'a constant for sufficient decrease';
    while 'alpha is not good enough'
        alpha = 'a shorter step length';
    end
    alpha_k = 'an alpha that is good enough';
end

% Function returning the value of the Rosenbrock function at x
function fval = f(x) % x is a vector of two variables
    fval = 'the Rosenbrock function';
end

% Function returning the value of the gradient at x
function grad = gradient(x)
    grad = [ 'df/dx1' ;
             'df/dx2' ];
end