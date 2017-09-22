close all
clear all

addpath(genpath(fullfile(pwd,'../common')))
addpath('common')
addpath('compute_nc')
addpath('compute')
addpath('compute_optimal')
addpath('compute/wavelet')

Sigmaids=8;%0:9;
Sigma_values = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024];
Nids=0;
C_true = loadbin('data/image_usi/solution.bin')';
width=1024;
height=512;
imagename='usi';

image_solution = reshape(C_true,width,height)';

%compute_nc_snr
%compute_nc_fourier
%compute_nc_fourierl2
%compute_nc_fouriersobolev
compute_nc_fouriervtr
%compute_nc_gauss
%compute_nc_median
%compute_nc_wiener
%compute_nc_denC2D
%compute_nc_denR2D
%compute_nc_denS2D

