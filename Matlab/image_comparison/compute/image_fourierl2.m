function [ Fout ] = image_fourierl2(F, s, lambda)

[height, width] = size(F);

% Kernel
x = [0:width/2-1, -width/2:-1];
y = [0:height/2-1, -height/2:-1];
[X,Y] = meshgrid(x,y);
h = exp( (-X.^2 - Y.^2)/(2*s^2) );
h = h/sum(h(:));

hF = real(fft2(h));
yF = fft2(F);

% We use this short hand for the filtering. Note that this is a symmetric operator.
Fout = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );

end

