function [ A,b,B,c,lb,ub ] = get_QP(C,Theta)

% in this implementation I assume 1D data
T = length(C);
K = length(Theta);

% create Hessian matrix
Ablock = sparse(2*diag(ones(T,1)) - diag(ones(T-1,1),1) - diag(ones(T-1,1),-1));
Ablock(1,1) = 1;
Ablock(T,T) = 1;
Ablock = sparse(Ablock);
A = sparse(K*T,K*T);
for k = 1:K
   A((k-1)*T+1:k*T,(k-1)*T+1:k*T) = Ablock; 
end

% matrix of equality constraints
B = sparse(T,K*T);
for k = 1:K
   B(:,(k-1)*T+1:k*T) = eye(T);
end
c = ones(T,1);

% lower bound
lb = zeros(K*T,1);
ub = ones(K*T,1);

% create RHS
b = zeros(K*T,1);
for t = 1:T
    for k = 1:K
        b((k-1)*T+t) = (C(t) - Theta(k))^2; % local model error
    end
end

% scale data
%A = (1/T)*A;
%b = (1/T)*b;

end

