close all
clear all

addpath('solvers')
addpath('functions')

global Mom    % yes, I am lazy

k = 2;
eps = 0;

% compute moments
%Mom = get_Mom(k);

% define some funny data
T = 10000;
mu = 2;
sigma = 3;
x_orig = mu + sqrt(sigma)*randn(T,1);

% x to [-1,1]
minx = min(x_orig);
maxx = max(x_orig);
ax = 2/(maxx-minx);
bx = 1 - ax*maxx;
x = ax.*x_orig + bx;

% compute moments
Mom = zeros(k,1);
for i = 1:k
    Mom(i) = sum(x.^i)/T;
    disp(['Mom_' num2str(i) ' = ' num2str(Mom(i))])
end


% some initial whatever
LM0 = zeros(k,1);

[LM,it] = mynewton(LM0,1e-6);
%[LM,it] = bb(LM0,1e-6);
%[LM,it] = sbb(LM0,1e-6);

% print solution
%LM'

% print C++ solution
for i=1:k
   disp(['Matlab_solution(' num2str(i-1) ') = ' num2str(LM(i)) ';']) 
end


