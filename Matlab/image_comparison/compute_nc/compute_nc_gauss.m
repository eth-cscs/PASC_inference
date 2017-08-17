disp('- computing gauss')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        N=[5 7 11];
        Sigma=[2 3 4];
        [F_gauss, Ferr_gauss, Nbest_gauss, Sigmabest_gauss] = optimal_gauss(image_data,image_solution, N, Sigma);
        disp(['   - best for N     = ' num2str(Nbest_gauss)])
        disp(['   - best for Sigma = ' num2str(Sigmabest_gauss)])
        disp(['   - error          = ' num2str(Ferr_gauss)])
    
        Ferr(n_mc_idx,i) = Ferr_gauss;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_gauss_' imagename '.mat'],'output');
