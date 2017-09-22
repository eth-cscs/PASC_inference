function [ output_args ] = check_constraints( gamma, K )

T = length(gamma)/K;
B = kron(ones(1,K),eye(T));

norm(B*gamma - ones(T,1))

end

