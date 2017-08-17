function e = double_dual1(s,x,t)

% EXAMPLE:
% s1 = double(imread('peppers.jpg'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% t = 0:5:45;
% e = double_dual1(s,x,t);
% plot(t,e);
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

N = length(t);
for k = 1:N
    y = doubledual_S1D(x,t(k));
    e(k) = sqrt(mean(mean((y-s).^2)));
end