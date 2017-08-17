disp('- computing fourier')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        S=[5 20 30 40 60 80];
        [F_fourier, Ferr_fourier, Sbest_fourier] = optimal_fourier(image_data,image_solution, S);
        disp(['     - best for S = ' num2str(Sbest_fourier)])
        disp(['     - error      = ' num2str(Ferr_fourier)])
    
        Ferr(n_mc_idx,i) = Ferr_fourier;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_fourier_' imagename '.mat'],'output');
