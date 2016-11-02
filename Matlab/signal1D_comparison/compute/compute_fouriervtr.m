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
    
    %%  TVR regularisation
    %  Through all noises
    for i=1:size(C,2)
        
        % f0 is a signal with noise which has to be filtered
        f0 = C(:,i)';
        if mod(size(f0,2),2) == 1
            f0 = f0(:,1:end-1);
        end
        
        % compute Fourier for different regularization parameters "s" and "lambda".
        lambda=[0.1 1000];
        S=[5 20];
        epsilon = [0.001 0.1];
        %figure;title(['i = ' num2str(i)]);
        clear F_err yF
        % for every combinations of parameters filter the signal
        for s=1:length(S)
            for l=1:length(lambda)
                for ee = 1:length(epsilon)
                    yF = compute_fourier_TVRBB_1D(f0, S(s), epsilon(ee), lambda(l), 50);

                    % compute the filtering errors F_err for different s
                    F_err(s,l,ee)=mean(abs(yF-C_true));
                end
            end
        end
        
        % compute minimal error
        [F_err_min(i,n_mc) s_best l_best ee_best] = minminmin(F_err);
        
        % provide information about best parameters
        disp(['  - best for Sigma=' num2str(Sigmaids(i)) ': s=' num2str(S(s_best)) ', lambda=' num2str(lambda(l_best)) ', epsilon=' num2str(epsilon(ee_best))])
        
        
    end
    
end

% prepare output and save it
output.F = F_err_min;

save(['results/output_fouriervtr_' nameextension '.mat'],'output');
