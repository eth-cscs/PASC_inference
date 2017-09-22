% compare denoising results of different methods on 100 signals
% this script does not compute, it only loads computed results
% to compute these resusts call functions in "compute" folder
%
% Lukas Pospisil, USI, Lugano 2017
%

close all
clear all

imagename = 'usi';
K = 2;

addpath(genpath(fullfile(pwd,'../common')))

% load computed data from Matlab
output_snr = load(['results_n100/output_snr_' imagename '.mat']);
output_fourier = load(['results_n100/output_fourier_' imagename '.mat']);
output_fourierl2 = load(['results_n100/output_fourierl2_' imagename '.mat']);
output_fouriersobolev = load(['results_n100/output_fouriersobolev_' imagename '.mat']);
output_fouriervtr = load(['results_n100/output_fouriervtr_' imagename '.mat']);
output_gauss = load(['results_n100/output_gauss_' imagename '.mat']);
output_median = load(['results_n100/output_median_' imagename '.mat']);
output_wiener = load(['results_n100/output_wiener_' imagename '.mat']);
output_denC2D = load(['results_n100/output_denC2D_' imagename '.mat']);
output_denR2D = load(['results_n100/output_denR2D_' imagename '.mat']);
output_denS2D = load(['results_n100/output_denS2D_' imagename '.mat']);

% load computed data form petsc
output_petsc = load_shortinfo_image(['results_n100/shortinfo_final.txt'],K);

% plot figure with errors
figure
hold on

title('Image Filtering Comparison');

xlabel('SNR');
ylabel('Filtering Error E_t[||x_{true}(t)-x_{filtered}(t)||]')
set(gca,'xscale','log','yscale','log')

x = output_snr.output.snr';
plot(x,output_fourier.output.Favg,'b.:','LineWidth',1.5,'MarkerSize',20);
plot(x,output_fourierl2.output.Favg,'r.:','LineWidth',1.5,'MarkerSize',20);
plot(x,output_fouriersobolev.output.Favg,'m.:','LineWidth',1.5,'MarkerSize',20);
plot(x,output_fouriervtr.output.Favg,'.:','Color',[0 0.4 0.0],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_gauss.output.Favg,'.:','Color',[0 0.8 0.4],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_median.output.Favg,'.:','Color',[0.8 0.0 0.4],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_wiener.output.Favg,'.:','Color',[0.8 0.4 0.0],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_denC2D.output.Favg,'.:','Color',[0.4 0.4 0.8],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_denR2D.output.Favg,'.:','Color',[0.4 0.8 0.8],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_denS2D.output.Favg,'.:','Color',[0.8 0.4 0.8],'LineWidth',1.5,'MarkerSize',20);
plot(x,output_petsc/(1000*600),'.:','Color',[0 0 0.4],'LineWidth',1.5,'MarkerSize',20);

%'Fourier VTR',
legend('Fourier', 'Fourier L2', 'Fourier Sobolev', 'Fourier TVR', 'Gauss', 'Median', 'Wiener', 'denC2D', 'denR2D', 'denS2D', 'H1FEM')

%axis([10^-4 10^0 10^-3 10^3])
hold off
