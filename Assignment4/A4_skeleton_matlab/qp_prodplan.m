% ********************************************************
% *                                                      *
% *      Optimalisering og regulering 					 *
% * 		?ving 3 Oppgave   V?r 2003					 *
% *                                                      *
% *      Bjarne Foss 1996                                *
% *                                                      *
% * qp_prodplan.m                                        *
% *                                                      *
% * m-file for calculating QP solution.                  *
% *                                                      *
% * Oppdated 10/1-2001 by Geir Stian Landsverk           *
% *                                                      *
% * Verified to work with MATLAB R2015a,                 *
% *    Andreas L. Fl?ten                                 *
% *                                                      *
% * Verified to work with MATLAB R2018b,                 *
% *    Joakim R. Andersen                                *
% *                                                      *
% ********************************************************


global XIT; % Storing iterations in a global variable
global IT;  % Storing number of iterations in a global variable
global x0;  % Used in qp1.m (line 216).

IT=1; XIT=[];

x0 = [0 0]'; % Initial value
vlb = [0 0]; % Lower bound on x
vub = [100 100];    % Upper bound on x

% min 0.5*x'*G*x + x'*c
%  x 
%
% s.t. A*x <= b 

% Objective
G = [0.8 0; 0 0.4];
c = [-3; -2];

% Ineq. constraints
A = -[-2 -1;
     -1 -3;
      1  0;
      0  1];
b = -[-8; -15; 0; 0];

options = optimset('LargeScale','Off');
x = quadprog1(G,c,A,b,[],[],vlb,vub,x0,options);

disp('Iteration sequence:')
disp(XIT');
disp('Solution:')
disp(x);
