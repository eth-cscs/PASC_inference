close all
clear all

%addpath(genpath(fullfile(pwd,'../common')))
addpath('util')

% number of clusters (for gamma0)
K = 2;
epssqr = 100;
Theta = [1,2];

% create original signal
C_small = [2*ones(1,60) ones(1,120) 2*ones(1,140) ones(1,230) 2*ones(1,150) ones(1,225) 2*ones(1,75)];
%C_small = [2*ones(1,6) ones(1,12) 2*ones(1,14) ones(1,23) 2*ones(1,15) ones(1,23) 2*ones(1,7)];
Sigma = 10^1;
nmb_of_repeat = 5;

% repeat small signal
C_true=repmat(C_small,1,nmb_of_repeat);
T = length(C_true);
disp([' - size of C_true: ' num2str(T)])

% signal with noise
C=C_true+sqrt(Sigma)*randn(size(C_true)); % add gauss noise
      
% generate gamma0
x0 = randn(K*T,1);

% create QP problem
[A,b,B,c,lb,ub] = get_QP(C,Theta);
A = epssqr*A;

%interior-point-convex
%x_ip = zeros(size(x0));
x_ip = test_ip(A,-b,B,c,lb,ub,x0);

%active-set
x_as = zeros(size(x0));
%x_as = test_as(A,-b,B,c,lb,ub,x0);


% reconstruct solution
C_reconstructed_ip = zeros(T,1);
C_reconstructed_as = zeros(T,1);
for k=1:K
   C_reconstructed_ip = C_reconstructed_ip + Theta(k)*x_ip((k-1)*T+1:k*T); 
   C_reconstructed_as = C_reconstructed_as + Theta(k)*x_as((k-1)*T+1:k*T); 
end

% plot problem
figure
hold on
plot(1:length(C),C,'b');
plot(1:length(C),C_true,'r');
plot(1:length(C),C_reconstructed_ip,'m');
%plot(1:length(C),C_reconstructed_as,'g');
%legend('with noise','exact', 'ip', 'as')
legend('with noise','exact', 'reconstructed')
xlabel('t')
ylabel('X(t)')
hold off


