function y = doubledualtree_i1D(w, J, Fsf, sf)

% Inverse 1-D Double-Density Dual-Tree Complex DWT
%
% USAGE:
%    y = dualtree_i1D(w, J, Fsf, sf)
% INPUT:
%    w - DWT coefficients
%    J - number of stages
%    Fsf - synthesis filters for the last stage
%    sf - synthesis filters for preceeding stages
% OUTUT:
%    y - output signal
% See also: doubledualtree_f1D
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

% Tree 1
y1 = w{J+1}{1};
for j = J:-1:2
   y1 = sfb3(y1, w{j}{1}, sf{1});
end
y1 = sfb3(y1, w{1}{1}, Fsf{1});

% Tree 2
y2 = w{J+1}{2};
for j = J:-1:2
   y2 = sfb3(y2, w{j}{2}, sf{2});
end
y2 = sfb3(y2, w{1}{2}, Fsf{2});

% normalization
y = (y1 + y2)/sqrt(2);