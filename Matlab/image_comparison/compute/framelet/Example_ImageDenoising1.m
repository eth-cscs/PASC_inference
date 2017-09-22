s1 = double(imread('peppers.jpg'));	% load image as a double
s = s1(:,:,3);                      % convert to a 2-D image
figure(1)                           % figure 1 displays the original image
imagesc(s)
colormap(gray)
axis image
title('Original Image')
x = s + 20*randn(size(s));          % add Gaussian noise to the original image
figure(2)                           % figure 2 displays the noisy image
imagesc(x)
colormap(gray)
axis image
title('Noisy Image')
T = 20;                             % choose a threshold of 20
y = double_S2D(x,T);             % denoise the noisy image using Double-Density Dual-Tree DWT
figure(3)                           % figure 3 diplays the denoised image
imagesc(y)
colormap(gray)
axis image
title('Denoised Image')