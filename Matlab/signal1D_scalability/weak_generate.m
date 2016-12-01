clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% number of clusters (for gamma0)
K = 2;

% number of nodes (for weak scalability)
N = [1,2,3,4];

% First we create the artificial time series signal C_true
C_small=[2*ones(1,60) ones(1,40) 2*ones(1,200) ones(1,50) 2*ones(1,150) ones(1,50) 2*ones(1,150) ones(1,225) 2*ones(1,75)];
Sigma=0.1*2.^20;
noise_small = sqrt(Sigma)*randn(size(C_small));
gamma0_small = randn(1,K*length(C_small));

weak_repeat_and_save('data/weak_10e4',C_small,noise_small,gamma0_small,10,N);
weak_repeat_and_save('data/weak_10e5',C_small,noise_small,gamma0_small,100,N);
weak_repeat_and_save('data/weak_10e6',C_small,noise_small,gamma0_small,1000,N);
weak_repeat_and_save('data/weak_10e7',C_small,noise_small,gamma0_small,10000,N);


