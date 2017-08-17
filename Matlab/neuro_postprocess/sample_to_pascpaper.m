% results for PASC paper
%
% Lukas Pospisil, USI, Lugano 2016
%


clear all
close all

% the name of input file with results
filename = 'test_edf_epssqr0.1';

addpath(genpath(fullfile(pwd,'../common')))
addpath('function')

filename_orig = ['results/' filename '_original.bin'];
filename_rec = ['results/' filename '_recovered.bin'];

% Load brain mesh
load('data/brain_coordinates_simple.mat','fv');

% Load 3D coordinates of electrodes
load('data/cap_coordinates_3D.mat','cap_coordinates_3D');

% Load data of time-series
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);

R = size(cap_coordinates_3D,1); % = 64, number of electrodes

% add some data to electrodes, each timestep stored in column
data_orig = reshape(vec_orig,R,length(vec_orig)/R);
data_rec = reshape(vec_rec,R,length(vec_rec)/R);
T = size(data_orig,2);

% plot some electrodes as 1D signal
plot_ids = [1,2,3];
if true
    figure
    
    for i = 1:length(plot_ids)
        subplot(length(plot_ids),1,i)
        hold on

        plot(1:T,data_orig(plot_ids(i),:),'b')
        plot(1:T,data_rec(plot_ids(i),:),'r')
        
        xlabel('$t$','Interpreter','latex')
        ylabel(['$X_{' num2str(plot_ids(i)) '}(t)$'],'Interpreter','latex')

        if i==1
            legend('original','recovered')
        end
        hold off
    end
end


% plot 3D brain
if false
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

    set(gcf,'Renderer','OpenGL');
    set(gcf,'RendererMode','manual');

    for t=10:10:1000
        % plot electrodes
        if false
            plot3(cap_coordinates_3D(:,1),cap_coordinates_3D(:,2),cap_coordinates_3D(:,3),'k*');
        end

        % map electrode data to brain vertex data
        brain_data = map_electrode_to_brain( fv.vertices, cap_coordinates_3D, data_rec(:,t) );

        disp(['t = ' num2str(t)])
    
        % update brain coloring
        update_brain(brain_patch, brain_data);

        % set color axis
        caxis([-100 100])
    
        pause
    end
end

hold off

    