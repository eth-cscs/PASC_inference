function strong_repeat_and_save(filename_begin,C_small,nmb_of_repeat,Sigma,K)

% First we create the artificial time series signal C_true
C_true=repmat(C_small,1,nmb_of_repeat);
disp([' - size of C_true: ' num2str(length(C_true))])

% save solution
filename = strcat(filename_begin, '_solution.bin');
savebin( filename, C_true );

% signal with noise
C=C_true+sqrt(Sigma)*randn(size(C_true)); % add gauss noise
       
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
hold off

end