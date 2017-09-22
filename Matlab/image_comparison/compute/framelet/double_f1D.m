function w = double_f1D(x, J, af)

% Forward Double-Density Discrete 1-D Wavelet Transform
%
% USAGE:
%    w = double_f1D(x, J, af)
% INPUT:
%    x - N-point vector, where
%            1) N is divisible by 2^J
%            2) N >= 2^(J-1)*length(af)
%    J - number of stages
%    af - analysis filters
%    af(:, 1) - lowpass filter (even length)
%    af(:, 2) - first highpass filter (even length)
%    af(:, 3) - second highpass filter (even length)
% OUTPUT:
%    w{j}, j = 1...J+1 - double-density DWT coefficients
% EXAMPLE:
%    [af, sf] = filters1;
%    x = rand(1,64);
%    w = double_f1D(x,3,af);
%    y = double_i1D(w,3,sf);
%    err = x - y; 
%    max(abs(err))
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

for k = 1:J
    [x w{k}] = afb3(x, af);
end

w{J+1} = x;