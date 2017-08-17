function [ x, it ] = bb( x0, myeps )

% set initial approximation
x_old = x0;
g_old = get_grad(x_old);

% first step
alpha = 1.0; % stepsize

x = x0 - alpha*g_old;
g = get_grad(x);

it = 1;

while true % stopping criteria using "break"
    s = x - x_old;
    y = g - g_old;
%    alpha = dot(s,y)/dot(s,s);
    alpha = dot(s,s)/dot(s,y);

    x_old = x;
    g_old = g;
    
    x = x - alpha*g;
    g = get_grad(x);
    
    it = it+1;

    disp(['it=' num2str(it) ': norm(g) = ' num2str(norm(g)) ', norm(s) = ' num2str(norm(s)) ', f(LM)=' num2str(get_f(x)) ]) 
    
    % stopping criteria
    if it > 1000 || norm(g) < myeps
       break; 
    end
end

% and that's all, folks!

end