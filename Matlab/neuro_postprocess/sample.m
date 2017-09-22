% plot data from electrodes to (into,onto) brain
%
% Lukas Pospisil, USI, Lugano 2016
%


clear all
close all

addpath('function')

% Load brain mesh
load('data/brain_coordinates_simple.mat','fv');

% Load 3D coordinates of electrodes
load('data/cap_coordinates_3D.mat','cap_coordinates_3D');

nmb_electrodes = size(cap_coordinates_3D,1); % = 64

% add some data to electrodes
electrode_data = zeros(nmb_electrodes,1);
for i = 1:nmb_electrodes
   electrode_data(i) = cap_coordinates_3D(i,1); % = x coordinate of electrode (test)
end

% map electrode data to brain vertex data
brain_data = map_electrode_to_brain( fv.vertices, cap_coordinates_3D, electrode_data );

figure
hold on

plot_brain(fv.vertices,fv.faces,brain_data);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');

colorbar

axis off
view([-135 35]);

set(gcf,'Renderer','OpenGL');
set(gcf,'RendererMode','manual');

hold off
