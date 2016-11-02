% generate and store random initial gamma0
%
% Lukas Pospisil, USI, Lugano 2016
%


clear all
close all

addpath('myfunctions')

K = 2; % number of clusters
n = 10000; % length of signal

gamma0 = randn(1,K*n); % TODO: project to feasible set?
savebin( 'data/signal1D_gamma0.bin', gamma0 );

