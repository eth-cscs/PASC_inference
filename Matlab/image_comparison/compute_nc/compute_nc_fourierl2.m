disp('- computing fourier L2')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        Lambda=[0.00001 0.0001 0.001 0.01 0.1 1 10];
        S=[3 5 20 30 40 60 80];
        [F_fourierl2, Ferr_fourierl2, Sbest_fourierl2, Lambdabest_fourierl2] = optimal_fourierl2(image_data,image_solution, S, Lambda);
        disp(['     - best for S      = ' num2str(Sbest_fourierl2)])
        disp(['     - best for Lambda = ' num2str(Lambdabest_fourierl2)])
        disp(['     - error           = ' num2str(Ferr_fourierl2)])
    
        Ferr(n_mc_idx,i) = Ferr_fourierl2;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_fourierl2_' imagename '.mat'],'output');
