%clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

%gamma = loadbin('results/entropy_small_epssqr1e+08_gamma.bin');
gamma = loadbin('entropy_small_gamma.bin');

%gamma = z_end;

K = 2;
T = length(gamma)/K;

disp(['T = ' num2str(T)])
disp(['K = ' num2str(K)])

figure
for k=1:K
    subplot(K,1,k)
    hold on
    plot(1:T,gamma((k-1)*T+1:1:k*T),'r')

    hold off
end


