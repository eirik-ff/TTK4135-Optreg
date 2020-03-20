function [x_opt, fval_opt, x_iter, f_iter, alpha] = NewtonsDescent(x0, f, G, H, maxiter, grad_tol)
% Inputs
% x0: initial point
% f: function handle
% G: function handle for the gradient of f
% H: function handle for the hessian of f
% maxiter: maximum number of iterations (default 10000)
% grad_tol: minimum norm value of gradient (default 1e-4)
% 
% Returns 
% x_opt: x^*
% fval_opt: f(x^*)
% x_iter: all iterates k -- column k contains x_k
% f_iter: a vector of all function values f(x_k)
% alpha: a vector of all step lenghts alpha_k

% Termination criteria
maxiter_default  = 10000;
grad_tol_default = 1e-4;

if nargin == 4
    maxiter = maxiter_default;
    grad_tol = grad_tol_default;
elseif nargin == 5
    grad_tol = grad_tol_default;
end

nx = size(x0,1); % Number of variables

% Declare some variables
x       = NaN(nx,maxiter);
p       = NaN(nx,maxiter);
grad    = NaN(nx,maxiter);
hessian = NaN(nx,nx,maxiter);
alpha   = NaN(1,maxiter);
fval    = NaN(1,maxiter);


% Do some calcualtions before the while loop. I.e., do the first iteration:
k = 1; % iteration number
x(:,k) = x0;
% look at the different functions below and finish them before continuing.
fval(k) = f(x0);
grad(:,k) = G(x0);
hessian(:,:,k) = H(x0);
p(:,k) = nd(grad(:,k), hessian(:,k)); 

alpha_0 = 1;

alpha(k) = linesearch(f, x(:,k), p(:,k), fval(k), grad(:,k), alpha_0);
x(:,k+1) = x(:,k) + alpha(k) * p(:,k);
grad(:,k+1) = G(x(:,k+1));
hessian(:,:,k+1) = H(x(:,k+1));
k = k + 1;

while k < maxiter && norm(grad(:,k),2) >= grad_tol
    fval(k) = f(x(:,k));    % Evaluate the Rosenbrock function
    p(:,k) = nd(grad(:,k), hessian(:,:,k));
    alpha_0 = 1;
    alpha(k) = linesearch(f, x(:,k), p(:,k), fval(k), grad(:,k), alpha_0); % Determine alpha using Alg. 3.1
    x(:,k+1) = x(:,k) + alpha(k) * p(:,k);
    grad(:,k+1) = G(x(:,k+1));
    hessian(:,:,k+1) = H(x(:,k+1));
    k = k + 1;
end
fval(k) = f(x(:,k)); % Final function value

% Delete unused space
x = x(:,1:k);
p = p(:,1:k);
grad = grad(:,1:k);
hessian = hessian(:,:,1:k);
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
function p = nd(grad, hessian)
    p = -hessian \ grad;
end

% Function implementing Algorithm 3.1, page 37 in N&W
function alpha_k = linesearch(f, xk, pk, fk, gradk, alpha_0)
    alpha = alpha_0;
    rho = 0.9;
    c1 = 0.5; % a constant for sufficient decrease
    while f(xk + alpha*pk) > fk + c1*alpha*(gradk.')*pk
        alpha = rho * alpha;
    end
    alpha_k = alpha;
end
