function y = double_S1D(x,T)

% x: noise signal
% T: threshold
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

[af, sf] = filters1;
J = 4;
w = double_f1D(x,J,af);
% loop through scales
for j = 1:J
    % loop through subbands
    for s = 1:2
       w{j}{s} = soft(w{j}{s},T);
    end
end
y = double_i1D(w,J,sf);