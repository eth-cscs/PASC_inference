function [ y0 ] = compute_fourier_Sobolev_1D(f0, s,lambda)

% get size of signal 
nx = length(f0);

% Kernel
x = [0:nx/2-1, -nx/2:-1];
h = exp( (-x.^2)/(2*s^2) );
h = h/sum(h(:));
hF = real(fft2(h));
yF = fft2(f0);

% Compute the Sobolev prior penalty |S| (rescale to [0,1]).
S = (x.^2 + x.^2)*(2/nx)^2;

% Perform the inversion.
y0 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );



end

