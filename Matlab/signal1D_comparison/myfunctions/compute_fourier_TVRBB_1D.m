function [ y0, E, it ] = compute_fourier_TVRBB_1D(f0, s, epsilon, lambda, niter)

% get size of signal 
nx = length(f0);

% Kernel
x = [0:nx/2-1, -nx/2:-1];
h = exp( (-x.^2)/(2*s^2) );
h = h/sum(h(:));
h = h';

% Initialization.
fTV = f0';

% Compute the TV norm, usefull to keep track of its decay through iterations.
%tv = sum(d(:));

% define iteration mapping
Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));

tau = 1.9 / ( 1 + lambda * 8 / epsilon);
alpha = tau;

Gr = mygrad(fTV);
d = sqrt( epsilon^2 + sum(Gr.^2,3) );
G = -mydiv( Gr./d );
e = Phi(fTV,h)-f0';
dd = Phi(e,h) + lambda*G;

% we use E for postprocessing ( to see descent )
clear E
it = 1;
for it = 1:niter
    % step
    fTVold = fTV;
    fTV = fTV - alpha*dd;
    
    % gradient
    ddold = dd;
    Gr = mygrad(fTV);
    d = sqrt( epsilon^2 + sum(Gr.^2,3) );
    G = -mydiv( Gr./d );
    e = Phi(fTV,h)-f0';
    dd = Phi(e,h) + lambda*G;
    
    s = fTV - fTVold;
    alpha = dot(s,s)/dot(s,dd-ddold);
    
    % energy
    E(it) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
%    disp([num2str(it) ',' num2str(E(it)) ',' num2str(norm(s)) ',' num2str(norm(dd-ddold)) ',' num2str(norm(dd))])

    if it>1
%        disp([num2str(it) ': ' num2str(norm(E(it)-E(it-1)))])
        if norm(E(it)-E(it-1)) < sqrt(length(fTV))
            break
        end
    end
end

% set output
y0 = fTV';

end

function g = mygrad(f)
    g = f([2:end end],:,:)-f;
end

function fx = mydiv(Px)
    fx = Px-Px([1 1:end-1],:,:);         
    fx(1,:,:)   = Px(1,:,:);        % boundary
    fx(end,:,:) = -Px(end-1,:,:);        
end
