function [ Fout ] = image_median(F, N)
% see https://ch.mathworks.com/help/images/noise-removal.html

Fout = medfilt2(F,[N N]);

end

