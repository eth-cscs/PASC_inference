% compare denoising results of different methods on 100 signals
% this script does not compute, it only loads computed results
% to compute these resusts call functions in "compute" folder
%
% (Fourier, FourierL2, FourierSobolev, FourierVTR, HMM, FEM-H1)
%
% Lukas Pospisil, USI, Lugano 2016
%

close all

addpath('myfunctions')

%nameextension = 'gauss';
nameextension = 'uniform';

% load computed data from Matlab
output_snr = load(['results_n100/output_snr_' nameextension '.mat']);
output_fourier = load(['results_n100/output_fourier_' nameextension '.mat']);
output_fourierl2 = load(['results_n100/output_fourierl2_' nameextension '.mat']);
output_fouriersobolev = load(['results_n100/output_fouriersobolev_' nameextension '.mat']);
output_fouriervtr = load(['results_n100/output_fouriervtr_' nameextension '.mat']);
output_hmm = load(['results_n100/output_hmm_' nameextension '.mat']);

% load computed data form petsc
output_petsc = load_shortinfo(['results_n100/shortinfo_final_' nameextension '.txt']);

% plot figure with errors
figure
hold on

title('Image Filtering Comparison');

xlabel('SNR');
ylabel('Filtering Error E_t[||x_{true}(t)-x_{filtered}(t)||]')
set(gca,'xscale','log','yscale','log')

x = mean(output_snr.output.snr');
plot(x,mean(output_fourier.output.F'),'b.:','LineWidth',1.5,'MarkerSize',20);
plot(x,mean(output_fourierl2.output.F'),'r.:','LineWidth',1.5,'MarkerSize',20);
plot(x,mean(output_fouriersobolev.output.F'),'m.:','LineWidth',1.5,'MarkerSize',20);
plot(x,mean(output_fouriervtr.output.F'),'.:','Color',[0 0.4 0.0],'LineWidth',1.5,'MarkerSize',20);
plot(x,mean(output_hmm.output.F'),'.:','Color',[0 0.8 0.4],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_petsc,'.:','Color',[0 0 0.4],'LineWidth',1.5,'MarkerSize',20);

legend('Fourier', 'Fourier L2', 'Fourier Sobolev', 'Fourier TVR', 'HMM', 'H1FEM')

axis([10^-4 10^0 10^-3 10^3])
hold off
