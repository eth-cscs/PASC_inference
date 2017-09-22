disp('- computing fourier Sobolev')

Ferr = zeros(length(Nids),length(Sigmaids));

for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);
    disp(['  - n_mc = ' num2str(n_mc) ])

    % for every sigma load C with noise
    for i = 1:length(Sigmaids)
        C_data = loadbin(['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'])';        
        image_data = reshape(C_data,width,height)';

        disp(['    - sigma = ' num2str(Sigmaids(i))])
        
        Lambda=[0.001 0.01 0.1 1 10 50 100];
        S=[3 5 20 30 40 60 80];
        [F_fouriersobolev, Ferr_fouriersobolev, Sbest_fouriersobolev, Lambdabest_fouriersobolev] = optimal_fouriersobolev(image_data,image_solution, S, Lambda);
        disp(['     - best for S      = ' num2str(Sbest_fouriersobolev)])
        disp(['     - best for Lambda = ' num2str(Lambdabest_fouriersobolev)])
        disp(['     - error           = ' num2str(Ferr_fouriersobolev)])
    
        Ferr(n_mc_idx,i) = Ferr_fouriersobolev;
    end
    
end

% prepare output and save it
output.Favg = mean(Ferr,1)';
output.Fmax = max(Ferr,1)';
output.Fmin = max(Ferr,1)';

save(['results_n100/output_fouriersobolev_' imagename '.mat'],'output');
