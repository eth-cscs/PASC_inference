clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% definition of this simple sample
n = 100;
K = 1;
mu = 0;
sigma = 1;
filename_begin = 'data/entropy_small';

C_true = mu + sigma*randn(n,1);

disp([' - size of C: ' num2str(length(C_true))])

% save solution
filename = strcat(filename_begin, '_solution.bin');
savebin( filename, C_true );

% signal with noise
C=C_true;%+sqrt(Sigma)*randn(size(C_true)); % add gauss noise
       
filename = strcat(filename_begin, '_data.bin');
savebin( filename, C );

% generate gamma0
gamma0 = randn(1,K*length(C));
filename = strcat(filename_begin, '_gamma0.bin');
savebin( filename, gamma0 );

figure
hold on
plot(1:length(C),C,'b');
plot(1:length(C),C_true,'r');
legend('with noise','exact')
xlabel('t')
ylabel('X(t)')
hold off
