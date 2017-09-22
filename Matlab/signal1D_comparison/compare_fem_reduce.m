% compare results of different FEM_REDUCE on 100 signals
% this script does not compute, it only loads computed results
% to compute these resusts call functions in "compute" folder
%
% Lukas Pospisil, USI, Lugano 2017
%

close all
clear all

addpath(genpath(fullfile(pwd,'../common')))

nameextension = 'gauss';
nameextension2 = 'gauss_1';
%nameextension = 'uniform';

% load original data (without noise) to get t
C_true = loadbin(['data/' nameextension '/signal1D_solution.bin']);
T = length(C_true);

% load computed data from Matlab
output_snr = load(['results/output_snr_' nameextension '.mat']);
output = load(['results_fem/' nameextension2 '.mat']);

%for i=1:length(output.fem_reduce)
%    output{i} = output_temp.output.F;
%end

% plot figure with errors
figure
hold on

%title('Image Filtering Comparison');

xlabel('signal noise ratio', 'Interpreter', 'latex', 'FontSize', 12);
%ylabel('Filtering Error E_t[||x_{true}(t)-x_{filtered}(t)||]', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('filtering error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca,'xscale','log','yscale','log')

x = mean(output_snr.output.snr');

mycolors{1} = [0 0 1];
mycolors{2} = [1 0 0];
mycolors{3} = [1 0 1];
mycolors{4} = [0 0.4 0];
mycolors{5} = [0 0.8 0.4];
mycolors{6} = [0 0 0.4];

for i=1:length(output.output.fem_reduce)
    plot(x,output.output.F_fem{i}/T,'.-','Color',mycolors{i},'LineWidth',2.5,'MarkerSize',20);
end

%h = legend('Fourier', 'Fourier L2', 'Fourier Sobolev', 'Fourier TVR', 'HMM', 'FEM-H1');
%set(h,'Interpreter','latex');
%set(h, 'FontSize', 12);

%axis([10^-4 10^0 10^-3 10^0])
hold off
