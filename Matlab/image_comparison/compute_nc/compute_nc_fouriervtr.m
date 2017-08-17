disp('- computing fourier VTR')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        Lambda=[100 1000];
        S=[10 20 50 80];
        Epsilon = [0.01 0.1 0.001];
        [F_fouriertvr, Ferr_fouriertvr, Sbest_fouriertvr, Lambdabest_fouriertvr, Epsilonbest_fouriertvr] = optimal_fouriertvr(image_data,image_solution, S, Lambda, Epsilon);
        disp(['   - best for S       = ' num2str(Sbest_fouriertvr)])
        disp(['   - best for Lambda  = ' num2str(Lambdabest_fouriertvr)])
        disp(['   - best for Epsilon = ' num2str(Epsilonbest_fouriertvr)])
        disp(['   - error            = ' num2str(Ferr_fouriertvr)])
    
        Ferr(n_mc_idx,i) = Ferr_fouriertvr;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_fouriervtr_' imagename '.mat'],'output');
