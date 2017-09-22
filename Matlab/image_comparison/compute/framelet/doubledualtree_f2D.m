function w = doubledualtree_f2D(x, J, Faf, af)

% Forward Real 2-D Double-Density Dual-Tree DWT
%
% USAGE:
%   w = doubledualtree_f2D(x, J, Faf, af)
% INPUT:
%   x - M by N array
%   J - number of stages
%   Faf - analysis filters for first stage
%   af - analysis filters for remaining stages
% OUPUT:
%   w{j}{d1}{d2} - DWT coefficients
%       j = 1..J, k = 1..2, d = 1..3
%   w{J+1}{k} - lowpass coefficients
%       k = 1..2
% % EXAMPLE:
%   x = rand(256);
%   J = 3;
%   [Faf, Fsf] = FSdoubledualfilt;
%   [af, sf] = doubledualfilt;
%   w = doubledualtree_f2D(x, J, Faf, af);
%   y = doubledualtree_i2D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

% normalization
x = x/sqrt(2);

% Tree 1
[x1 w{1}{1}] = afb3_2D(x, Faf{1});      % stage 1
for j = 2:J
    [x1 w{j}{1}] = afb3_2D(x1, af{1});  % remaining stages
end
w{J+1}{1} = x1;                       % lowpass subband

% Tree 2
[x2 w{1}{2}] = afb3_2D(x, Faf{2});      % stage 1
for j = 2:J
    [x2 w{j}{2}] = afb3_2D(x2, af{2});  % remaining stages
end
w{J+1}{2} = x2;                       % lowpass subband

% sum and difference
for j = 1:J
    for m = 1:8
        A = w{j}{1}{m};
        B = w{j}{2}{m};
        w{j}{1}{m} = (A+B)/sqrt(2);
        w{j}{2}{m} = (A-B)/sqrt(2);
    end
end