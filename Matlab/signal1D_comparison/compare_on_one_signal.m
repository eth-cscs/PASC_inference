% compare denoising results of different methods on one particular signal
% (Fourier, FourierL2, FourierSobolev, FourierVTR, HMM, FEM-H1)
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all
close all

addpath(genpath(fullfile(pwd,'../common')))

% load HMM toolbox
addpath('../HMM/HMMall/HMM')
addpath('../HMM/HMMall/KPMstats')
addpath('../HMM/HMMall/KPMtools')
addpath('../HMM/HMMall/netlab3.3')

% compute signal with "uniform" or "gauss" noise?
nameextension = 'uniform';
%nameextension = 'gauss';

Sigmaid=5; % id of sigma
n_mc=1; % id of generated noise
% i.e. load signal 'data/'+nameextension+'/signal1D_id'+n_mc+'_idSigma'+Sigmaid+'.bin'

% load original data (without noise)
C_true = loadbin(['data/' nameextension '/signal1D_solution.bin']);

% load data with noise
C_noise = loadbin(['data/' nameextension '/signal1D_id' num2str(n_mc) '_idSigma' num2str(Sigmaid) '.bin']);

%% compute Fourier
disp('- compute fourier')
f0 = C_noise';
if mod(size(f0,2),2) == 1
 f0 = f0(:,1:end-1);
end
% compute fourier for different sizes of window "s".
S=[5 20 30 40 60 80];
F_err_best = Inf;
C_fourier = f0;
for s=1:length(S)
    yF = compute_fourier_1D(f0, S(s));
    F_err = mean(abs(yF-C_true'));
    if F_err < F_err_best
        F_err_best = F_err_best;
        s_best = s;
        C_fourier = yF;
    end
end
disp(['  - best: s=' num2str(S(s_best))])

%% compute FourierL2
disp('- compute fourier L2')
lambda=[0.1 1 10 50 100];
S=[5 20 30 40 60 80];
F_err_best = Inf;
C_fourierl1 = f0;
for s=1:length(S)
	for l=1:length(lambda)
        yF = compute_fourier_L2_1D(f0, S(s),lambda(l));
        F_err = mean(abs(yF-C_true'));
        if F_err < F_err_best
            F_err_best = F_err;
            s_best = s;
            l_best = l;
            C_fourierl2 = yF;
        end
    end
end
disp(['  - best: s=' num2str(S(s_best))  ', lambda=' num2str(lambda(l_best))])

%% compute FourierSobolev
disp('- compute fourier sobolev')
% compute Fourier for different regularization parameters "s" and "lambda".
lambda=[0.1 1 10 50 100];
S=[5 20 30 40 60 80];
F_err_best = Inf;
C_fouriersobolev = f0;
for s=1:length(S)
    for l=1:length(lambda)
        yF = compute_fourier_Sobolev_1D(f0, S(s),lambda(l));
        F_err = mean(abs(yF-C_true'));
        if F_err < F_err_best
            F_err_best = F_err;
            s_best = s;
            l_best = l;
            C_fouriersobolev = yF;
        end
    end
end
disp(['  - best: s=' num2str(S(s_best))  ', lambda=' num2str(lambda(l_best))])

%% compute FourierVTR
disp('- compute fourier VTR')
lambda=[0.1 1000];
S=[5 20];
epsilon = [0.001 0.1];
F_err_best = Inf;
C_fouriervtr = f0;
for s=1:length(S)
	for l=1:length(lambda)
        for ee = 1:length(epsilon)
            yF = compute_fourier_TVRBB_1D(f0, S(s), epsilon(ee), lambda(l), 50);
            F_err = mean(abs(yF-C_true'));
            if F_err < F_err_best
                F_err_best = F_err;
                s_best = s;
                l_best = l;
                ee_best = ee;
                C_fouriervtr = yF;
            end
        end
    end
end
disp(['  - best: s=' num2str(S(s_best))  ', lambda=' num2str(lambda(l_best)) ', epsilon=' num2str(epsilon(ee_best))])

%% compute HMM
disp('- compute HMM')
%% if HMM was already computed, then use it, otherwise compute it
result_hmm_mat = ['results/hmm_' nameextension '/signal1D_' nameextension '_id' num2str(n_mc) '_idSigma_' num2str(Sigmaid) '_recovered.mat'];
if exist(result_hmm_mat, 'file') == 2
    C_hmm_load = load(result_hmm_mat);
    C_hmm = C_hmm_load.C_hmm;
    disp(['  - loaded from ' result_hmm_mat])
else    
    F_err_best = Inf;
    C_hmm = f0;
    for annealing_id=1:10
        disp([' - computing annealing: ' num2str(annealing_id)])
        % the algorithm uses rand to assemble initial guess, 
        % therefore annealing means one random realisation
        yF = compute_hmm_1D(f0);
        F_err = mean(abs(yF-C_true'));
        if F_err < F_err_best
            F_err_best = F_err;
            annealing_id_best = annealing_id;
            C_hmm = yF;
        end
    end
    disp(['  - best : annealing_id=' num2str(annealing_id_best)])
    save(result_hmm_mat,'C_hmm');
end


%% load h1fem solution computed with our library
disp('- load FEM-H1 solution')
C_h1fem = loadbin(['results/' nameextension '/signal1D_' nameextension '_id' num2str(n_mc) '_idSigma_' num2str(Sigmaid) '_recovered.bin']);
  

%% show solution
figure
title([nameextension ', id = ' num2str(n_mc) ', idSigma = ' num2str(Sigmaid)])
hold on
plot(1:length(C_noise),C_noise,'-','Color',[0.7 0.7 0.7],'LineWidth',1);
plot(1:length(C_true),C_true,'-','Color',[0.2 0.2 0.2],'LineWidth',1.5);
plot(1:length(C_fourier),C_fourier,'b-','LineWidth',1);
plot(1:length(C_fourierl2),C_fourierl2,'r-','LineWidth',1);
plot(1:length(C_fouriersobolev),C_fouriersobolev,'m-','LineWidth',1);
plot(1:length(C_fouriervtr),C_fouriervtr,'-','Color',[0 0.4 0.0],'LineWidth',1);
plot(1:length(C_hmm),C_hmm,'-','Color',[0 0.8 0.4],'LineWidth',1);
plot(1:length(C_h1fem),C_h1fem,'-','Color',[0 0 0.4],'LineWidth',1);
legend('with noise', 'without noise', 'Fourier', 'Fourier L2', 'Fourier Sobolev', 'Fourier TVR', 'HMM', 'FEM-H2')

axis([1 length(C_true) -2 5])
hold off
