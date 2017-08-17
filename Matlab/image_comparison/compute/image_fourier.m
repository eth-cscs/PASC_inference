function [ Fout ] = image_fourier(F, s)

[height, width] = size(F);

% Kernel
x = [0:width/2-1, -width/2:-1];
y = [0:height/2-1, -height/2:-1];
[X,Y] = meshgrid(x,y);
h = exp( (-X.^2 - Y.^2)/(2*s^2) );
h = h/sum(h(:));

% We use this short hand for the filtering. Note that this is a symmetric operator.
Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));

% Apply the filter.
Fout = Phi(F,h);

end

