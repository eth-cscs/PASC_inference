disp('- computing median')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        N=[5 10 20 30];
        [F_median, Ferr_median, Nbest_median] = optimal_median(image_data,image_solution, N);
        disp(['   - best for N     = ' num2str(Nbest_median)])
        disp(['   - error          = ' num2str(Ferr_median)])
    
        Ferr(n_mc_idx,i) = Ferr_median;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_median_' imagename '.mat'],'output');
