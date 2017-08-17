function [ x, it ] = mynewton( x0, myeps )

% set initial approximation
x = x0;
it = 0;

alpha = 1.0; % magic dumping coefficient (0,1]

while true % stopping criteria using "break"
    H = get_hess( x );
    g = get_grad( x );
    
    eigH = eig(H);
    disp(['condition number: ' num2str(max(eigH)/min(eig(H))) ])
    disp(['min             : ' num2str(min(eig(H))) ])
    disp(['max             : ' num2str(max(eig(H))) ])
    
    eigH'
    
    deltax = -H\g; % H=H', H >= 0 => CG method?   
    
    x = x + alpha*deltax;

    it = it+1;

    disp(['it=' num2str(it) ': norm(deltax) = ' num2str(norm(deltax)) ', f(LM)=' num2str(get_f(x)) ]) 

    % stopping criteria
    if it > 100 || norm(deltax) < myeps
       break; 
    end
end

% and that's all, folks!

end