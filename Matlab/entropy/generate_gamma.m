clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% generate gamma0
T = 1000;
gamma0 = [ones(1,500) 0*ones(1,T-500) 0*ones(1,500) ones(1,T-500) ];
filename = 'big_test_gamma0.bin';
savebin( filename, gamma0 );

