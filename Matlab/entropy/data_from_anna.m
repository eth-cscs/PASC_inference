clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% Anna: sample from paper
load('data/xt.mat');
y = xt;

% Anna: toy sample (K=3, there should be switch in t=1000, t=2000 )
%load('data/y.mat');
filename_begin = 'data/xt';

figure
hold on
plot(1:length(y),y,'b')
hold off

filename = strcat(filename_begin, '_data.bin');
savebin( filename, y );


for k=1:5
    gamma0 = randn(1,k*length(y));
    filename = strcat(filename_begin, ['_K' num2str(k) '_gamma0.bin']);
    savebin( filename, gamma0 );
end




