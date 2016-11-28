clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% number of clusters (for gamma0)
K = 2;

% First we create the artificial time series signal C_true
C_small=[2*ones(1,60) ones(1,40) 2*ones(1,200) ones(1,50) 2*ones(1,150) ones(1,50) 2*ones(1,150) ones(1,225) 2*ones(1,75)];
Sigma=0.1*2.^20;

strong_repeat_and_save('data/strong_10e4',C_small,10,Sigma,K);
strong_repeat_and_save('data/strong_10e5',C_small,100,Sigma,K);
strong_repeat_and_save('data/strong_10e6',C_small,1000,Sigma,K);
strong_repeat_and_save('data/strong_10e7',C_small,10000,Sigma,K);


