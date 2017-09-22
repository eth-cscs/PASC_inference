function [ gamma ] = get_random_gamma0( T, K )
%GET_RANDOM_GAMMA0
%

% generate initial random feasible gamma
gamma = rand(K*T,1); % in (0,1)
gamma_sum = zeros(T,1);
for k=1:K
    gamma_sum = gamma_sum+gamma((k-1)*T+1:k*T);
end
for k=1:K
    gamma((k-1)*T+1:k*T) = gamma((k-1)*T+1:k*T)./gamma_sum; % to obtain sum_{k=1:K} gamma_k = 1
end


end

