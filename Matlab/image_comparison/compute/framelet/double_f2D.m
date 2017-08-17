function w = double_f2D(x, J, af1, af2)

% Forward Double-Density Discrete 2-D Wavelet Transform
%
% USAGE:
%   w = double_f2D(x, J, af)
% INPUT:
%   x - N by M matrix
%      1) M, N are both even
%      2) min(M,N) >= 2^(J-1)*length(af)
%   J - number of stages
%  af - analysis filters
% OUPUT:
%   w - cell array of wavelet coefficients
% EXAMPLE:
%   [af, sf] = filters1;
%   x = rand(128,64);
%   J = 3;
%   w = double_f2D(x,J,af);
%   y = double_i2D(w,J,sf);
%   err = x - y; 
%   max(max(abs(err)))
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

if nargin < 4
    af2 = af1;
end

for k = 1:J
    [x w{k}] = afb3_2D(x, af1, af2);
end

w{J+1} = x;