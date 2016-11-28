function [ y0 ] = compute_fourier_L2_1D(f0, s,lambda)

% get size of image 
nx = length(f0);

% Kernel
x = [0:nx/2-1, -nx/2:-1];
h = exp( (-x.^2)/(2*s^2) );
h = h/sum(h(:));
hF = real(fft2(h));
yF = fft2(f0);
% We use this short hand for the filtering. Note that this is a symmetric operator.
y0 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
% Apply the filter.

end

