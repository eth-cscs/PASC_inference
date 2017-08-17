function y = double_S2D(x,T)

% x: noise signal
% T: threshold 
% EXAMPLE
%  s1 = double(imread('st.tif'));
%  s = s1(:,:,3);
%  x = s + 20*randn(size(s));
%  T = 35;
%  y = double_denS2D(x,T);
%  imagesc(y)
%  colormap(gray)
%  axis image
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

[af, sf] = filters1;
J = 3; % number of stages
w = double_f2D(x,J,af);
% loop thru scales
for j = 1:J
    % loop thru subbands
    for s = 1:8
        w{j}{s} = soft(w{j}{s},T);
    end
end
y = double_i2D(w,J,sf); 