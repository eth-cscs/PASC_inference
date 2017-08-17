%clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

nameextension = 'gauss_1';
%nameextension = 'uniform';



% load shortinfo file
filename = ['results_fem/' nameextension '.txt' ];

%1col ?fem_reduce,idfile,idSigma,K,epssqr,
%6col abserr,energy,Theta0,Theta1,it all, 
%11col t all, t gamma update, t gamma solve, t theta update, t theta solve, 
%16col SPGQP it, SPGQP hessmult, SPGQP t all, SPGQP t project, SPGQP t matmult, 
%21col SPGQP t dot, SPGQP t update, SPGQP t stepsize, SPGQP t fs, SPGQP t other, 
%26col SPGQP fx, SPGQP fx_linear, SPGQP fx_quadratic, SPGQP_sum it, SPGQP_sum hessmult, 
%31col SPGQP_sum t all, SPGQP_sum t project, SPGQP_sum t matmult, SPGQP_sum t dot, SPGQP_sum t update, 
%36col SPGQP_sum t stepsize, SPGQP_sum t fs, SPGQP_sum t other,        

M = csvread(filename,1,0);

output = [];
output.fem_reduce = sort(unique(M(:,1)));
output.file_id = sort(unique(M(:,2)));
output.sigma_id = sort(unique(M(:,3)));

output.F_fem = [];
for femid = 1:length(output.fem_reduce)
    disp(['- femid=' num2str(femid)])

    errs = zeros(1,length(output.sigma_id));
    for sigmaid = 1:length(output.sigma_id)
        disp([' - sigmaid=' num2str(sigmaid)])
    
        abserr = zeros(1,length(output.file_id));

        for n_mc=1:length(output.file_id)
            aa = M(M(:,1)==output.fem_reduce(femid) & M(:,2)==output.file_id(n_mc) & M(:,3)==output.sigma_id(sigmaid),6);
            if length(aa) > 0
                abserr(n_mc) = min(aa);
            else
                disp(['missing: signal1D_id' num2str(output.file_id(n_mc)) '_idSigma' num2str(output.sigma_id(sigmaid)) '_fem' num2str(output.fem_reduce(femid))])
            end
        end
        
        errs(sigmaid) = sum(abserr)/length(output.file_id);
    end

    output.F_fem{femid} = errs;
end

% save output
save(['results_fem/' nameextension '.mat'],'output');

