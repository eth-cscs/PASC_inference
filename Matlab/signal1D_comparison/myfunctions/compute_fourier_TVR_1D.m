function [ y0, E ] = compute_fourier_TVR_1D(f0, s, epsilon, lambda, niter)

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

% we use E for postprocessing ( to see descent )
clear E
for i=1:niter
    % Compute the gradient of the smoothed TV functional.
    Gr = mygrad(fTV);
    d = sqrt( epsilon^2 + sum(Gr.^2,3) );

%    G = -div( Gr./repmat(d, [1 1 2])  );
    G = -mydiv( Gr./d );

    % step
    e = Phi(fTV,h)-f0';
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
    % energy
    E(i) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
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
