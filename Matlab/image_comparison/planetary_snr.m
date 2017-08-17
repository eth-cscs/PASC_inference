clear all

addpath('common')
addpath('compute')
addpath('compute/wavelet')
addpath('compute_optimal')

C_true = loadbin(['data/planetary_solution.bin'])';
C_data = loadbin(['data/planetary_data.bin'])';

[snrvalue, snrvalue_sigma] = snr(C_data, C_true, 0.064);
disp(['SNR = std(solution)/std(data) = ' num2str(snrvalue) ])
disp(['SNR = std(solution)/sigma     = ' num2str(snrvalue_sigma) ])

