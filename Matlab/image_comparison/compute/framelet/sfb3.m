function y = sfb3(lo, hi, sf)

% Synthesis Filter Bank 
% (consisting of three filters)
%
% USAGE:
%    y = sfb3(lo, hi, sf)
% INPUT:
%     lo - low frqeuency input
%     hi - cell array containing high frequency inputs
%       1) hi{1} - first high frequency input
%       2) hi{2} - second high frequency input
%     sf - synthesis filters
%    sf(:, 1) - lowpass filter (even length)
%    sf(:, 2) - first highpass filter (even length)
%    sf(:, 3) - second highpass filter (even length)
% OUTPUT:
%    y - output signal
% See also: afb3
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

N = 2*length(lo);
L = length(sf);
lo  = upfirdn(lo, sf(:,1), 2, 1);
hi{1} = upfirdn(hi{1}, sf(:,2), 2, 1);
hi{2} = upfirdn(hi{2}, sf(:,3), 2, 1);
y = lo + hi{1} + hi{2};
y(1:L-2) = y(1:L-2) + y(N+[1:L-2]);
y = y(1:N);
y = cshift(y, 1-L/2);