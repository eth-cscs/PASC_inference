function y = doubledualtree_i2D(w, J, Fsf, sf)

% Inverse Real 2-D Double-Density Dual-Tree DWT
% 
% USAGE:
%   y = doubledualtree_i2D(w, J, Fsf, sf)
% INPUT:
%   J - number of stages
%   Fsf - synthesis filters for first stage
%   sf -  synthesis filters for remaining stages
% OUPUT:
%   y - output array
% See also: doubledualtree_f2D
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

% sum and difference
for k = 1:J
    for m = 1:8
        A = w{k}{1}{m};
        B = w{k}{2}{m};
        w{k}{1}{m} = (A+B)/sqrt(2);
        w{k}{2}{m} = (A-B)/sqrt(2);
    end
end

% Tree 1
y1 = w{J+1}{1};
for j = J:-1:2
   y1 = sfb3_2D(y1, w{j}{1}, sf{1});
end
y1 = sfb3_2D(y1, w{1}{1}, Fsf{1});

% Tree 2
y2 = w{J+1}{2};
for j = J:-1:2
   y2 = sfb3_2D(y2, w{j}{2}, sf{2});
end
y2 = sfb3_2D(y2, w{1}{2}, Fsf{2});

% normalization
y = (y1 + y2)/sqrt(2);