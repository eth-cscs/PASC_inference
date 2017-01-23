function [ x fx it ] = sqpqp(A,b,x0, K, myeps)
% implementation of SPG-QP, solve 1/2*x'*A*x - b'*x st. Bx=c, x>=l
% where B,c,l are from FEM-H1 (based on K)
%
% INPUT:
% A - SPS matrix
% b - RHS vector
% x0 - initial approximation
% K - number of clusters (defines equality constraints)
% myeps - accuracy
%
% OUTPUT:
% x - solution
% it - number of main iterations
%

% projection of initial approximation
x = get_projection(x0,K);
  
% SPG-QP settings (magic parameters)
m = 20;
gamma = 0.9; 
sigma2 = 1.0;
maxit = 100; % don't worry, we improve solution in next outer loop

it = 0;

% 0. Initialization
alpha_bar = 2.0;

g = A*x - b;

f_old = Inf;
f = get_f( x, g, b);

fs = f*ones(m,1);
alpha_bb = alpha_bar;

while it < maxit
    d = get_projection(x-alpha_bb*g,K) - x;
    dd = dot(d,d);
    
    Ad = A*d;
    dAd = dot(Ad,d);
    f_max = max(fs);
    
    xi = (f_max - f)/dAd;
    beta_bar = -dot(g,d)/dAd;
    beta_hat = gamma*beta_bar + sqrt(gamma^2*beta_bar^2 + 2*xi);
    
    beta = min([sigma2,beta_hat]);

    x = x + beta*d;
    g = g + beta*Ad;

    f_old = f;
    f = get_f( x, g, b);
    fs(1:end-1) = fs(2:end);
    fs(end) = f;
    
    alpha_bb = dd/dAd;

    it = it + 1;
    
    if abs(f-f_old) < myeps*norm(b)
        break;
    end

end

% return last function value
fx = f;

end

% compute function value 
function fx = get_fbrute(x,A,b)
    fx = 0.5*dot(A*x,x) - dot(x,b);
end

% compute function value from gradient (without matrix-vector multiplication)
function [ fx ] = get_f( x, g, b)
    fx = 1/2*dot(g-b,x);
end

% project point to feasible set
function Px = get_projection(x,K)
   % find Px = arg min || y - x || s.t. y >= 0 and Bx = c 
   T = length(x)/K;
   
   Px = x;
   for i = 1:T
       range = i:T:i+(K-1)*T;
       
       xlocal = x(range);
       Pxlocal = get_projection_local(xlocal); 
       Px(range) = Pxlocal;
       
   end
end

% project point onto simplex
function [Px] = get_projection_local(x)
    n = length(x); % = K

    [y,idx] = sort(x); 
       
	t_hat = 0;
       
	i = n-1;
    while i >= 1 
        ti = (sum(y(i+1:n))-1)/(n-i);
        if ti >= y(i)
            t_hat = ti;
            break
        else
            i = i - 1;
        end
    end
       
	if i == 0
        t_hat = (sum(y(1:n))-1)/n;
    end
       
    Px = max(x-t_hat,zeros(n,1));
end

