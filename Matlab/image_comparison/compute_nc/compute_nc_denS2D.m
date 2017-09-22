disp('- computing S2D')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        T=[0.01 0.05 0.1 0.5 1 2 5 10 15 20 25 30];
        [F_denS2D, Ferr_denS2D, Tbest_denS2D] = optimal_denS2D(image_data,image_solution, T);
        disp(['   - best for T     = ' num2str(Tbest_denS2D)])
        disp(['   - error          = ' num2str(Ferr_denS2D)])
    
        Ferr(n_mc_idx,i) = Ferr_denS2D;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_denS2D_' imagename '.mat'],'output');
