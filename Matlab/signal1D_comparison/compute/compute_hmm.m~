clear all
close all

addpath('myfunctions')
addpath('../HMM/HMMall/HMM')
addpath('../HMM/HMMall/KPMstats')
addpath('../HMM/HMMall/KPMtools')
addpath('../HMM/HMMall/netlab3.3')

% this computation takes quite a long time, I am interested in it to
% estimate computing time of larger problems
tic

%nameextension = 'gauss';
nameextension = 'uniform';

fileid = '';
Sigmaids=1:10;
N_mc=1:100;

C_true = loadbin(['data/' nameextension '/signal1D_solution.bin'])';

for n_mcit=1:length(N_mc)
    n_mc = N_mc(n_mcit);
    n_mc

    % for every sigma prepare C with noise
    for i = 1:length(Sigmaids)
       C(:,i) = loadbin(['data/' nameextension '/signal1D_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
    end
    
    %%  HMM
    %  Through all noises
    for i=1:size(C,2)
        
        % f0 is a signal with noise which has to be filtered
        f0 = C(:,i)';

        clear F_err
        for annealing_id=1:10
            disp([' - computing annealing: ' num2str(annealing_id)])
            % the algorithm uses rand to assemble initial guess, 
            % therefore annealing means one random realisation
            yF = compute_hmm_1D(f0);
            F_err(annealing_id)=mean(abs(yF-C_true));
        end
        
        % compute the minimal filtering error F_err_min
        [F_err_min(i,n_mc) annealing_id_best] = min(F_err);

        % provide information about best parameters
        disp(['  - best for Sigma_' num2str(Sigmaids(i)) ': annealing_id=' num2str(annealing_id_best)])
        
    end
    
end

% prepare output and save it
output.F = F_err_min;
output.N_mc = N_mc;
output.Sigmaids = Sigmaids;

save(['results/output_hmm_' nameextension fileid '.mat'],'output');

disp(['computed in ' num2str(toc) 's'])