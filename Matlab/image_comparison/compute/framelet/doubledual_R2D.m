function y = doubledual_R2D(x,T);

% % Example
% s1 = double(imread('peppers.jpg'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% T = 10;
% y = doubleden_R2D(x,T);
% imagesc(y)
% colormap(gray)
% axis image
% sqrt(mean(mean((y-s).^2)))
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

[Faf, Fsf] = FSdoubledualfilt;
[af, sf] = doubledualfilt;
J = 3;
w = doubledualtree_f2D(x,J,Faf,af);
% loop thru scales:
for j = 1:J
    % loop thru subbands
    for s1 = 1:2
        for s2 = 1:8
            w{j}{s1}{s2} = soft(w{j}{s1}{s2},T);
        end
    end
end
y = doubledualtree_i2D(w,J,Fsf,sf);