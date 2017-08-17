%
% Lukas Pospisil, USI, Lugano 2017
%

close all
clear all

imagename = 'usi';

addpath(genpath(fullfile(pwd,'../common')))

% load computed data form petsc
[output_petsc0, fem_reduces] = load_shortinfo_image_fem(['results_fem/shortinfo_final0.txt']);
[output_petsc1, fem_reduces] = load_shortinfo_image_fem(['results_fem/shortinfo_final1.txt']);

width=1024;
height=512;
%size_reduced=width*height*fem_reduces.^2;

% plot figure with errors
figure
hold on

title('FEM reduction Comparison');

xlabel('fem reduction');
ylabel('Filtering Error E_t[||x_{true}(t)-x_{filtered}(t)||]')
set(gca,'yscale','log')

plot(fem_reduces,output_petsc0/(1000*600),'.:','Color',[0 0 0.4],'LineWidth',1.5,'MarkerSize',20);
plot(fem_reduces,output_petsc1/(1000*600),'.:','Color',[0.4 0 0],'LineWidth',1.5,'MarkerSize',20);

legend('FEM-SUM', 'FEM-HAT')

%axis([0 1 10^-3 10^3])
hold off
