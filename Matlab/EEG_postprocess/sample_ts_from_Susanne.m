% plot ts data from electrodes to (into,onto) brain
%
% Lukas Pospisil, USI, Lugano 2016
%


clear all
close all

addpath('function')

% Load brain mesh
load('data/brain_coordinates.mat','fv');

% Load 3D coordinates of electrodes
load('data/cap_coordinates_3D.mat','cap_coordinates_3D');

% Load data of time-series
load('data/input_for_clustering_problem_from_Susanne.mat','x');

R = size(cap_coordinates_3D,1); % = 64, number of electrodes

% add some data to electrodes, each timestep stored in column
electrode_data = reshape(x,R,length(x)/R);


figure
hold on

% initial plot, returns pointer to brain object
brain_patch = plot_brain(fv.vertices,fv.faces,zeros(size(fv.vertices)));

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');

colorbar

axis off
view([-135 35]);


for t=1:10
    % map electrode data to brain vertex data
    brain_data = map_electrode_to_brain( fv.vertices, cap_coordinates_3D, electrode_data(:,10) );

    disp(['t = ' num2str(t)])
    
    % update brain coloring
    update_brain(brain_patch, brain_data);

    pause
end

hold off

    