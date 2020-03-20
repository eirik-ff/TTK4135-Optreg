function [xopt, fval, kiters, xiters, Witers, piters, alphaiters lambdaiters] = activesetConvexQPSolver(G, c, A, b, x0, W0, zeroTol)
% Implementation of Algorithm 16.3 on page 472 of Nocedal & Wright's book
% "Numerical Optimization". 
%
% NOTE: NOT THOROUGHLY TESTED, MIGHT NOT GIVE CORRECT ANSWER. IT'S
% RECOMMENDED TO USE MATLABS quadprog INSTEAD. THIS IMPLEMENTATION IS JUST
% FOR FUN. 
%
% This is an Active-Set method for primal convex QP problems of the form
%
%   min 1/2 * x' * G * x + c' * x
%    x
% 
%   s.t. A*x >= b
%        Aeq*x = beq  (NOT SUPPORTED)
%
% The algorithm uses x0 as the initial point. It must be feasible. 
% W0 is the initial working set. Might be empty. 
% zeroTol: used to determine how close pk must be to zero to accept it as
%          zero
%
% Return variables are:
%   xopt: optimal point
%   fval: objective function value in xopt
%   kiters: the k iterations (k = 0,1,2,...)
%   xiters: the intermediate x-values for each iteration
%   Witers: the working set Wk for each iteration. Entries are NaN if
%      working set is empty.
%   piters: the step direction for each iteration
%   alphaiters: the alpha step length for each iteration (if pk = 0, this
%      value doesn't make sense)
%   lambdaiters: the lagrange multiplier lambda for each iteration. Entries
%      are NaN if the working set is empty. 

nearZero = @(x, tol) abs(x) < tol ;

kiters = [];
xiters = [];
Witers = [];
piters = [];
alphaiters = [];
lambdaiters = [];

problemSolved = 0;  
xopt = zeros(size(x0));

k = 0;
xk = x0;
xk1 = xk;
Wk = W0;
Wk1 = Wk; 
while not(problemSolved)
    xk = xk1;
    Wk = Wk1;
    
    gk = G * xk + c;
    K = [G, -A(Wk,:)'; A(Wk,:), zeros(size(A(Wk,:),1))];
    plambdak = K \ [-gk; zeros(length(Wk))];
    
    pk = plambdak(1:end-length(Wk));
    lambdak = plambdak(end-length(Wk)+1:end);
    
    alphaks = ones(size(b));
    for i = 1:size(b,1)
        if A(i,:) * pk < 0  % only check when ai^T pk < 0
            alphaks(i) = (b(i) - A(i,:) * xk) / (A(i,:) * pk);
        end  % if
    end  % for
    [alphak, i] = min(alphaks);
    
    
    if all(nearZero(pk, zeroTol))
        if all(lambdak >= 0)
            problemSolved = 1;
            xopt = xk;
        else 
            [~,j] = min(lambdak);
            j = Wk(j);
            Wk1 = Wk(Wk~=j);
            xk1 = xk;
        end  % if
    else  % pk != 0
        xk1 = xk + alphak * pk;
        
        if alphak < 1
            Wk1 = [Wk, i];
        else
            Wk1 = Wk;
        end  % if
    end  % if
    
    kiters = [kiters, k];
    xiters = [xiters, xk];
    piters = [piters, pk];
    alphaiters = [alphaiters, alphak];
    if isempty(Wk)
        Witers = [Witers, NaN];
        lambdaiters = [lambdaiters, NaN];
    else
        Witers = [Witers, Wk];
        lambdaiters = [lambdaiters, lambdak];
    end  % if    
    
    k = k + 1;
end  % while 

fval = 1/2 * xopt' * G * xopt + c' * xopt;
end  % function


