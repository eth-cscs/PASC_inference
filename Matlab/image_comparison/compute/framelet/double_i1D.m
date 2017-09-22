function y = double_i1D(w, J, sf)

% Inverse Double-Density Discrete 1-D Wavelet Transform
%
% USAGE:
%     y = double_i1D(w, J, sf)
% INPUT:
%     w - wavelet coefficients
%     J - number of stages
%    sf - synthesis filters
% OUTPUT:
%     y - output signal
%
% See also: double_f1D
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

y = w{J+1};

for k = J:-1:1
   y = sfb3(y, w{k}, sf);
end