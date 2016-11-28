function [ y0 ] = compute_hmm_1D(f0)

% the input has to be row
if size(f0,1) > size(f0,2)
   f0 = f0';
end

T = size(f0,2); % length of time-series
O = size(f0,1); % dimension of data
K = 2; % number of clusters
M = 1; % number of gauss mixtures
cov_type = 'full'; % type of used covariance matrix

% initial guess of parameters
prior0 = normalise(rand(K,1));
transmat0 = mk_stochastic(rand(K,K));
[mu0, Sigma0] = mixgauss_init(K*M, f0, cov_type);

% prepare initial guess
mu0 = reshape(mu0, [O K M]);
Sigma0 = reshape(Sigma0, [O O K M]);
mixmat0 = mk_stochastic(rand(K,M));

% run EM algorithm
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
        mhmm_em(f0, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 200);

% compute log-likelihood value
loglik = mhmm_logprob(f0, prior1, transmat1, mu1, Sigma1, mixmat1);

% deal with viterbi to obtain clustering functions
B = mixgauss_prob(f0, mu1, Sigma1, mixmat1);
[G] = viterbi_path(prior1, transmat1, B);

y0 = zeros(size(G));
y0(G == 1) = mu1(1);
y0(G == 2) = mu1(2);


end

