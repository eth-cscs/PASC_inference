function y = doubledual_C2D(x,T)

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

[Faf, Fsf] = FSdoubledualfilt;
[af, sf] = doubledualfilt;
I = sqrt(-1);
J = 3;
w = cplxdoubledual_f2D(x,J,Faf,af);
% loop thru scales
for j = 1:J
    % loop thru subbands
    for s1 = 1:2
        for s2 = 1:8
            C = w{j}{1}{s1}{s2} + I*w{j}{2}{s1}{s2};
            C = soft(C,T);
            w{j}{1}{s1}{s2} = real(C);
            w{j}{2}{s1}{s2} = imag(C);
        end
    end
end
y = cplxdoubledual_i2D(w,J,Fsf,sf); 