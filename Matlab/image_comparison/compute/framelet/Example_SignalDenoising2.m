s1 = double(imread('peppers.jpg'));	% load image as a double
s = s1(:,4,3);                      % convert 3-D image to a 1-D signal
subplot (3,1,1), plot(s)            % plot the original signal
title('Original Signal')
axis([0 512 -100 300])
x = s + 20*randn(size(s));          % add Gaussian noise to the original signal
subplot (3,1,2), plot(x)            % plot the noisy signal
title('Noisy Signal')
axis([0 512 -100 300])
T = 20;                             % choose a threshold of 20
y = doubledual_S1D(x,T);            % denoise the noisy image using Double-Density Dual-Tree DWT
subplot (3,1,3), plot(y)            % plot the denoised signal
title('Denoised Signal')
axis([0 512 -100 300])