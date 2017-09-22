function [ Fout ] = image_wiener(F, N)
% see https://ch.mathworks.com/help/images/noise-removal.html

Fout = wiener2(F,[N N]);

end

