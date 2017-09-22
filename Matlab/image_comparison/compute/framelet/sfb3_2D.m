function y = sfb3_2D(lo, hi, sf1, sf2)

% 2-D Synthesis Filter Bank
%
% USAGE:
%   y = sfb3_2D(lo, hi, sf1, sf2);
% INPUT:
%    lo - lowpass subband
%    hi - cell array containing eight highpass subbands
%   sf1 - synthesis filters for the columns
%   sf2 - synthesis filters for the rows
% OUTPUT:
%   y - output array
% See also: afb3_2D
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 


if nargin < 4
    sf2 = sf1;
end

% filter along rows
L = sfb3_2D_A(    lo, hi{1}, hi{2}, sf2, 2);
H1 = sfb3_2D_A(hi{3}, hi{4}, hi{5}, sf2, 2);
H2 = sfb3_2D_A(hi{6}, hi{7}, hi{8}, sf2, 2);

% filter along columns
y = sfb3_2D_A(L, H1, H2, sf1, 1);


% LOCAL FUNCTION

function y = sfb3_2D_A(L, H1, H2, sf, d)

% 2-D Synthesis Filter Bank
% (along one dimension only)
%
% USAGE:
% y = sfb2D_A(lo, hi1, hi2, sf, d);
%
% INPUT:
%    lo, hi1, hi2 - one lowpass and two highpass subbands
%    sf - synthesis filters
%     d - dimension of filtering
% OUTPUT:
%    y - output signal
%
% see also: afb3_2D_A


lpf  = sf(:, 1);    % lowpass filter
hpf1 = sf(:, 2);    % first highpass filter
hpf2 = sf(:, 3);    % second highpass filter

if d == 2
   L = L';
   H1 = H1';
   H2 = H2';
end

N = 2*size(L,1);
len = length(sf);
y = upfirdn(L, lpf, 2, 1) + upfirdn(H1, hpf1, 2, 1) + upfirdn(H2, hpf2, 2, 1);
y(1:len-2, :) = y(1:len-2, :) + y(N+[1:len-2], :);
y = y(1:N, :);
y = cshift2D(y, 1-len/2);

if d == 2
   y = y';
end
