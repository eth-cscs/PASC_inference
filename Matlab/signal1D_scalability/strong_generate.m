clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% number of clusters (for gamma0)
K = 2;

% First we create the artificial time series signal C_true
C_small=[2*ones(1,600) ones(1,1200) 2*ones(1,1400) ones(1,2300) 2*ones(1,1500) ones(1,2250) 2*ones(1,750)];
%C_small=[2*ones(1,60) ones(1,40) 2*ones(1,150) ones(1,100) 2*ones(1,150) ones(1,100) 2*ones(1,100) ones(1,225) 2*ones(1,75)];
Sigma=10^1;

%strong_repeat_and_save('data/strong_10e4',C_small,1,Sigma,K);
%strong_repeat_and_save('data/strong_10e5',C_small,10,Sigma,K);
%strong_repeat_and_save('data/strong_10e6',C_small,100,Sigma,K);
%strong_repeat_and_save('data/strong_10e7',C_small,1000,Sigma,K);
strong_repeat_and_save('data/strong_10e8',C_small,10000,Sigma,K);
%strong_repeat_and_save('data/small',C_small,1,Sigma,K);


