% create petsc file with lorem ipsum text image and add noise
%
% Lukas Pospisil, USI, Lugano 2017
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

noise_coeff = 0.2;

image_width = 1000;
image_height = 1000;

filename_in = 'data/image_lorem/loremipsum.png';
filename_image_out = 'output/image_lorem/loremipsum.bin';
filename_graph_out = 'output/image_lorem/graph.bin';

[image_vec, nodes_vec] = saveimage( filename_in, filename_image_out, filename_graph_out, noise_coeff );
    

