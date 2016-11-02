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
    
    %%  L2-regularized Fourier filtering
    %  Through all noises
    for i=1:size(C,2)
        
        % f0 is a signal with noise which has to be filtered
        f0 = C(:,i)';
        if mod(size(f0,2),2) == 1
            f0 = f0(:,1:end-1);
        end
        
        % compute Fourier for different regularization parameters "s" and "lambda".
        lambda=[0.1 1 10 50 100];
        S=[5 20 30 40 60 80];

        clear F_err yF
        % for every combinations of parameters filter the signal
        for s=1:length(S)
            for l=1:length(lambda)
                yF = compute_fourier_L2_1D(f0, S(s),lambda(l));

                % compute the filtering errors F_err for different s
                F_err(s,l)=mean(abs(yF-C_true));
            end
        end
        
        % compute the minimal filtering error F_err_min
        [F_err_min(i,n_mc) s_best l_best] = minmin(F_err);
        
        % provide information about best parameters
        disp(['  - best for Sigma=' num2str(Sigmaids(i)) ': s=' num2str(S(s_best)) ', lambda=' num2str(lambda(l_best))])
        
    end
    
end

% prepare output and save it
output.F = F_err_min;

save(['results/output_fourierl2_' nameextension '.mat'],'output');
