% plot the results - original signal, solution and recoved signal

% include functions
addpath(genpath(fullfile(pwd,'../common')));

filename_rec = 'results_fem/signal1D_id10_idSigma2_fem001_recovered.bin';
vec_rec = loadbin(filename_rec);

figure 
hold on
plot(1:length(vec_rec),vec_rec,'r')
hold off