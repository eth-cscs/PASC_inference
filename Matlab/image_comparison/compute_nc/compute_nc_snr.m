for n_mc_idx=1:length(Nids)
    n_mc = Nids(n_mc_idx);

    % for every sigma load C with noise and compute SNR
    for i = 1:length(Sigmaids)
       filename = ['data/image_' imagename '/' imagename '_id' num2str(n_mc) '_idSigma' num2str(Sigmaids(i)) '.bin'];
       
       disp(['processing: ' filename])
        
       C(:,i) = loadbin(filename)';        

       SNR(i,n_mc_idx)=std(C_true)/std(C(:,i));
       
       if n_mc_idx==1
           SNR_sigma(i)=std(C_true)/Sigma_values(i);
       end
    end
end

% prepare output and save it
output.snr = SNR;
output.snr_sigma = SNR_sigma;

save(['results_n100/output_snr_' imagename '.mat'],'output');
