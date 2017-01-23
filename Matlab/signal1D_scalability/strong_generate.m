clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% number of clusters (for gamma0)
K = 2;

% First we create the artificial time series signal C_true
C_small=[2*ones(1,600) ones(1,400) 2*ones(1,2000) ones(1,500) 2*ones(1,1500) ones(1,500) 2*ones(1,1500) ones(1,2250) 2*ones(1,750)];
Sigma=10^3;

%strong_repeat_and_save('data/strong_10e4',C_small,1,Sigma,K);
%strong_repeat_and_save('data/strong_10e5',C_small,10,Sigma,K);
%strong_repeat_and_save('data/strong_10e6',C_small,100,Sigma,K);
strong_repeat_and_save('data/strong_10e7',C_small,1000,Sigma,K);
%strong_repeat_and_save('data/strong_10e8',C_small,10000,Sigma,K);


