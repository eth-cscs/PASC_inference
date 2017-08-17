function [ Fout ] = image_gauss(F, N, sigma)
% see https://stackoverflow.com/questions/27499057/how-do-i-create-and-apply-a-gaussian-filter-in-matlab-without-using-fspecial-im

%disp(['  - N = ' num2str(N) ', sigma = ' num2str(sigma)])

% Generate Gaussian mask
ind = -floor(N/2) : floor(N/2);
[X Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = h / sum(h(:));

% Convert filter into a column vector
h = h(:);

% Filter our image
I_pad = padarray(F, [floor(N/2) floor(N/2)]);
C = im2col(I_pad, [N N], 'sliding');
C_filter = sum(bsxfun(@times, C, h), 1);
Fout = col2im(C_filter, [N N], size(I_pad), 'sliding');

end

