function re = real_den2(s,x,t)

% % Example
% s1 = double(imread('peppers.jpg'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% t = 0:5:45;
% e = rdoubleden2(s,x,t);
% plot(t,e);
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

N = length(t);
for k = 1:N
    y = doubledual_R2D(x,t(k));
    re(k) = sqrt(mean(mean((y-s).^2)));
end