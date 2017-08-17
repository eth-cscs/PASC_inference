function [ Fout ] = image_fouriersobolev(F, s, lambda)

[height, width] = size(F);

% Kernel
x = [0:width/2-1, -width/2:-1];
y = [0:height/2-1, -height/2:-1];
[X,Y] = meshgrid(x,y);
h = exp( (-X.^2 - Y.^2)/(2*s^2) );
h = h/sum(h(:));

hF = real(fft2(h));
yF = fft2(F);

% Compute the Sobolev prior penalty |S| (rescale to [0,1]).
S = (X.^2 + Y.^2)*(4/(width*height));

% Perform the inversion.
Fout = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );

end

