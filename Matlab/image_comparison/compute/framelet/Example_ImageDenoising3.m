s1 = double(imread('peppers.jpg'));	% load image as a double
s = s1(:,:,3);                      % convert to a 2-D image
figure(1)                           % display original image
imagesc(s)
colormap(gray)
axis image
title('Original Image')
x = s + 20*randn(size(s));          % add Gaussian noise to image
figure(2)                           % display noisy image
imagesc(x)
colormap(gray)
axis image
title('Noisy Image')
T = 15;                             % choose a threshold of 15
y = doubledual_C2D(x,T);            % denoise image using Double-Density Dual-Tree Complex DWT
figure(3)                           % diplay denoised image
imagesc(y)
colormap(gray)
axis image
title('Denoised Image')