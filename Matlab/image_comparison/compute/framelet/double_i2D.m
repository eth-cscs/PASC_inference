function y = double_i2D(w, J, sf1, sf2)

% Inverse 2-D Double-Density Discrete Wavelet Transform
%
% USAGE:
%   y = double_i2D(w, J, sf)
% INPUT:
%   w - wavelet coefficients
%   J - number of stages
%  sf - synthesis filters
% OUTPUT:
%   y - output array
% See also: double_f2D
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

if nargin < 4
    sf2 = sf1;
end

y = w{J+1};

for k = J:-1:1
   y = sfb3_2D(y, w{k}, sf1, sf2);
end