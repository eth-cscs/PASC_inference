disp('- computing wiener')

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
        [F_wiener, Ferr_wiener, Nbest_wiener] = optimal_wiener(image_data,image_solution, N);
        disp(['   - best for N     = ' num2str(Nbest_wiener)])
        disp(['   - error          = ' num2str(Ferr_wiener)])
    
        Ferr(n_mc_idx,i) = Ferr_wiener;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_wiener_' imagename '.mat'],'output');
