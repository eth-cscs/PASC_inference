function [ x, it ] = sbb( x0, myeps )

m = 1;

% set initial approximation
x_old = x0;
g_old = get_grad(x_old);

% first step
alpha = 1.0; % stepsize

x = x0 - alpha*g_old;
g = get_grad(x);

it = 1;

f = get_f(x);
fs = f*ones(m,1);

while true % stopping criteria using "break"
    s = x - x_old;
    y = g - g_old;
%    alpha = dot(s,y)/dot(s,s);
    alpha = dot(s,s)/dot(s,y);
    d = -alpha*g;
    [beta, it_gll] = gll(x,d,fs);

    x_old = x;
    g_old = g;
    x = x + beta*d;
    g = get_grad(x);
    f = get_f(x);

    fs(1:end-1) = fs(2:end);
    fs(end) = f;
    
    it = it+1;

    disp(['it=' num2str(it) ': norm(g)=' num2str(norm(g)) ', norm(s)=' num2str(norm(s)) ', it_gll=' num2str(it_gll) ', f(LM)=' num2str(f) ]) 
    
    % stopping criteria
    if it > 1000 || norm(g) < myeps
       break; 
    end
end

% and that's all, folks!

end


function [beta, it_gll] = gll(x_0,d,fs)
    gamma = 0.5; % magic parameter
    sigma_1 = 0;
    sigma_2 = 1 - sigma_1;
    
    it_gll = 0;
    
    f_max = max(fs);
    f_0 = get_f(x_0);
    g_0 = get_grad(x_0);
    
    beta = 1; % initial step size
    x_temp = x_0 + beta*d;
    f_temp = get_f(x_temp);

    delta = dot(g_0,d);
    
    while f_temp > f_max + gamma*beta*delta
        beta_temp = -0.5*beta^2*delta/(f_temp - f_0 - beta*delta);
        
        if beta_temp >= sigma_1 && beta_temp <= sigma_2*beta
            beta = beta_temp;
        else
            beta = beta/2;
        end
        
        x_temp = x_0 + beta*d;
        f_temp = get_f(x_temp);
        
        it_gll = it_gll + 1;
    end
end    
    