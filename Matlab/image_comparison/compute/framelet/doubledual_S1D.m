function y = doubledual_S1D(x,T)

% x - noise signal
% T - threshold
% 
% FRAMELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% 

[Faf, Fsf] = FSdoubledualfilt;
[af, sf] = doubledualfilt;
J = 4;
w = doubledualtree_f1D(x,J,Faf,af);
% loop thru scales:
for j = 1:J
    % loop thru subbands
    for s1 = 1:2
        for s2 = 1:2
            w{j}{s1}{s2} = soft(w{j}{s1}{s2},T);
        end
    end
end
y = doubledualtree_i1D(w,J,Fsf,sf);