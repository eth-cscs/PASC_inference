clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% number of clusters (for gamma0)
K = 2;

% First we create the artificial time series signal C_true
%C_small=[2*ones(1,600) ones(1,1200) 2*ones(1,1400) ones(1,2300) 2*ones(1,1500) ones(1,2250) 2*ones(1,750)];
C_small=[2*ones(1,60) ones(1,120) 2*ones(1,140) ones(1,230) 2*ones(1,150) ones(1,225) 2*ones(1,75)];
C_small=[2*ones(1,6) ones(1,12) 2*ones(1,14) ones(1,23) 2*ones(1,15) ones(1,22) 2*ones(1,8)];

%Sigma=10^1;
Sigma=10^-2;

strong_repeat_and_save('data/strong_10e2',C_small,1,Sigma,K);
%strong_repeat_and_save('data/strong_10e3',C_small,1,Sigma,K);
%strong_repeat_and_save('data/strong_10e4',C_small,1,Sigma,K);
%strong_repeat_and_save('data/strong_10e5',C_small,10,Sigma,K);
%strong_repeat_and_save('data/strong_10e6',C_small,100,Sigma,K);
%strong_repeat_and_save('data/strong_10e7',C_small,1000,Sigma,K);
%strong_repeat_and_save('data/strong_10e8',C_small,10000,Sigma,K);
%strong_repeat_and_save('data/small',C_small,1,Sigma,K);


