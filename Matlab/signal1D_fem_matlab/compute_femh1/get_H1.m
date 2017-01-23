function [ H, B ] = get_H1( T,K )
%GET_H1 

disp(' - creating matrices')

Hblock = sparse(2*diag(ones(T,1)) - diag(ones(T-1,1),1) - diag(ones(T-1,1),-1));
Hblock(1,1) = 1;
Hblock(end,end) = 1;
Hblock = sparse(Hblock);
H = sparse(K*T,K*T);
for k = 1:K
   H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = Hblock; 
end

B = sparse(T,K*T);
for k = 1:K
   B(:,(k-1)*T+1:k*T) = eye(T);
end


end

