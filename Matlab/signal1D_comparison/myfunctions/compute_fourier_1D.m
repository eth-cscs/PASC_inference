function [ y0 ] = compute_fourier_1D(f0, s)

% get size of image 
nx = length(f0);

% Kernel
x = [0:nx/2-1, -nx/2:-1];
h = exp( (-x.^2)/(2*s^2) );
h = h/sum(h(:));

% We use this short hand for the filtering. Note that this is a symmetric operator.
Phi = @(x,h)real(ifft(fft(x).*fft(h)));

% Apply the filter.
y0 = Phi(f0,h);

end

