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

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
       C(:,i) = loadbin(['data/' nameextension '/signal1D_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
    end
    
    %%  Fourier filtering
    %  Through all noises
    for i=1:size(C,2)
        
        % f0 is a signal with noise which has to be filtered
        f0 = C(:,i)';
        if mod(size(f0,2),2) == 1
            f0 = f0(:,1:end-1);
        end
        
        % compute fourier for different sizes of window "s".
        S=[5 20 30 40 60 80];

        clear F_err
        for s=1:length(S)
            yF = compute_fourier_1D(f0, S(s));
            F_err(s)=mean(abs(yF-C_true));
        end
        
        % compute the minimal filtering error F_err_min
        [F_err_min(i,n_mc) s_best] = min(F_err);
        
        % provide information about best parameters
        disp(['  - best for Sigma_' num2str(Sigmaids(i)) ': s=' num2str(S(s_best))])
    end
    
end

% prepare output and save it
output.F = F_err_min;

save(['results/output_fourier_' nameextension '.mat'],'output');
