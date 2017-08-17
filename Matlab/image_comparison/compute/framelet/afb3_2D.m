function [lo, hi] = afb3_2D(x, af1, af2)

% 2-D Analysis Filter Bank
%
% USAGE:
% [lo, hi] = afb3_2D(x, af1, af2);
% INPUT:
%    x - NxM matrix, where min(N,M) > 2*length(filter)
%           (N, M are even)
%    af1 - analysis filter for the columns
%    af2 - analysis filter for the rows
%    af(:, 1) - lowpass filter
%    af(:, 2) - first highpass filter
%    af(:, 3) - second highpass filter
% OUTPUT:
%     lo - lowpass subband
%     hi - cell array containing the eight following hipass subbandshi1, hi2 - one lowpass and two highpass subbands
%         1) hi{1} = 'lo hi1' subband
%         2) hi{2} = 'lo hi2' subband
%         3) hi{3} = 'hi1 lo' subband
%         4) hi{4} = 'hi1 hi1' subband
%         5) hi{5} = 'hi1 hi2' subband
%         6) hi{6} = 'hi2 lo' subband
%         7) hi{7} = 'hi2 hi1' subband
%         8) hi{8} = 'hi2 hi2' subband
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

if nargin < 3
   af2 = af1;
end

% filter along columns
[L, H1, H2] = afb3_2D_A(x, af1, 1);

% filter along rows
[lo,    hi{1}, hi{2}] = afb3_2D_A(L, af2, 2);
[hi{3}, hi{4}, hi{5}] = afb3_2D_A(H1, af2, 2);
[hi{6}, hi{7}, hi{8}] = afb3_2D_A(H2, af2, 2);


% LOCAL FUNCTION

function [L, H1, H2] = afb3_2D_A(x, af, d)

% 2-D Analysis Filter Bank
% (along one dimension only)
%
% USAGE:
% [lo, hi1, hi2] = afb3_2D_A(x, af, d);
% INPUT:
%    x - NxM matrix, where min(N,M) > 2*length(filter)
%           (N, M are even)
%    af - analysis filter for the columns
%    af(:, 1) - lowpass filter
%    af(:, 2) - first highpass filter
%    af(:, 3) - second highpass filter
%    d - dimension of filtering (d = 1 or 2)
% OUTPUT:
%     lo, hi1, hi2 - one lowpass and two highpass subbands
%
% EXAMPLE:
% x = rand(32,64);
% [af, sf] = filters1;
% [lo, hi1, hi2] = afb3_2D_A(x, af, 1);
% y = sfb2D_A(lo, hi1, hi2, sf, 1);
% err = x - y;
% max(max(abs(err)))

lpf  = af(:, 1);    % lowpass filter
hpf1 = af(:, 2);    % first highpass filter
hpf2 = af(:, 3);    % second highpass filter

if d == 2
   x = x';
end

N = size(x,1);
len = size(af,1)/2;
x = cshift2D(x,-len);

L = upfirdn(x, lpf, 1, 2);
L(1:len, :) = L(1:len, :) + L([1:len]+N/2, :);
L = L(1:N/2, :);

H1 = upfirdn(x, hpf1, 1, 2);
H1(1:len, :) = H1(1:len, :) + H1([1:len]+N/2, :);
H1 = H1(1:N/2, :);

H2 = upfirdn(x, hpf2, 1, 2);
H2(1:len, :) = H2(1:len, :) + H2([1:len]+N/2, :);
H2 = H2(1:N/2, :);

if d == 2
   L = L';
   H1 = H1';
   H2 = H2';
end