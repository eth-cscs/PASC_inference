% generate and store signal with noise in PETSc format
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

nameextension = 'gauss';
%nameextension = 'uniform';

% first we create the artificial time series signal C_true
C_true=[2*ones(1,60) ones(1,40) 2*ones(1,200) ones(1,50) 2*ones(1,150) ones(1,50) 2*ones(1,150) ones(1,225) 2*ones(1,75)];
C_true=repmat(C_true,1,10);
Sigma=0.1*2.^[1:2:20];

% save solution (signal without noise)
filename = ['data/' nameextension '/signal1D_solution.bin'];
savebin( filename, C_true );

% N_mc sets the number of Monte Carlo samples of the time series
N_mc=100;
for n_mc=1:N_mc
    n_mc

    % for every sigma prepare C with noise
    for i=1:length(Sigma)
        if strcmp(nameextension,'gauss')
            C(:,i)=C_true+sqrt(Sigma(i))*randn(size(C_true)); % gaussian
        end
        if strcmp(nameextension,'uniform')
            C(:,i)=C_true+sqrt(Sigma(i))*(rand(size(C_true)) - 0.5 ); % uniform
        end
        
        filename = ['data/' nameextension '/signal1D_id' num2str(n_mc) '_idSigma' num2str(i) '.bin'];
        savebin( filename, C(:,i) );
        
    end
    
end
