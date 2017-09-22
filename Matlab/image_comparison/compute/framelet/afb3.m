function [lo, hi] = afb3(x, af)

% Analysis Filter Bank 
% (consisting of three filters)
%
% USAGE:
%    [lo, hi] = afb3(x, af)
% INPUT:
%     x - N-point vector, where
%            1) N is even
%            2) N >= length(af)
%    af - analysis filters
%    af(:, 1) - lowpass filter (even length)
%    af(:, 2) - first highpass filter (even length)
%    af(:, 3) - second highpass filter (even length)
% OUTPUT:
%     lo - low frequecy output
%     hi - cell array containing high frequency outputs
%       1) hi{1} - first high frequency output
%       2) hi{2} - second high frequency output
% EXAMPLE:
%    [af, sf] = filters1;
%    x = rand(1,64);
%    [lo, hi] = afb3(x, af);
%    y = sfb(lo, hi, sf);
%    err = x - y; 
%    max(abs(err))
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

N = length(x);
L = length(af)/2;
x = cshift(x,-L);

% lowpass filter
lo = upfirdn(x, af(:,1), 1, 2);
lo(1:L) = lo(N/2+[1:L]) + lo(1:L);
lo = lo(1:N/2);

% first highpass filter
hi{1} = upfirdn(x, af(:,2), 1, 2);
hi{1}(1:L) = hi{1}(N/2+[1:L]) + hi{1}(1:L);
hi{1} = hi{1}(1:N/2);

% second highpass filter
hi{2} = upfirdn(x,af(:,3), 1, 2);
hi{2}(1:L) = hi{2}(N/2+[1:L]) + hi{2}(1:L);
hi{2} = hi{2}(1:N/2);