% create petsc file with "alzheimer" image data
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('common')

% input
C_noise_medium = load('data/image_alzheimer/C_noise_medium');

% output
filename_image_out_medium = 'output/image_alzheimer/C_noise_medium_w1683_h374.bin';
saveimage_frommat( C_noise_medium.C_noise, filename_image_out_medium);

