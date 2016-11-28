% create petsc file with USI text image and add noise
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

noise_coeffs = [0.0 0.2 0.4 1.0];
noise_coeffs_print = {'00' '02' '04' '10'};

usi_width = [250 500 1000];
usi_height = [150 300 600];

for i = 1:length(usi_width)
    for j = 1:length(noise_coeffs)
        disp(['exporting: '  num2str(usi_width(i)) 'x' num2str(usi_height(i))])
    
        name_part = ['usi_' num2str(usi_width(i)) '_' num2str(usi_height(i))];
        filename_in = ['data/image_usi/' name_part '.png'];
        filename_image_out = ['output/image_usi/' name_part '_' noise_coeffs_print{j} '.bin'];
        filename_graph_out = ['output/image_usi/graph_' name_part '.bin'];

        [image_vec, nodes_vec] = saveimage( filename_in, filename_image_out, filename_graph_out, noise_coeffs(j) );
    
    end
end

