clear all
close all

addpath('myfunctions')

%nameextension = 'gauss';
nameextension = 'uniform';

Sigmaids=1:10;
N_mc=100;

C_true = loadbin(['data/' nameextension '/signal1D_solution.bin'])';

for n_mc=1:N_mc
    n_mc

    % for every sigma load C with noise and compute SNR
    for i = 1:length(Sigmaids)
       C(:,i) = loadbin(['data/' nameextension '/signal1D_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
       SNR(i,n_mc)=1/(max(C(:,i))-min(C(:,i)));
    end
end

% prepare output and save it
output.snr = SNR;

save(['results/output_snr_' nameextension '.mat'],'output');
