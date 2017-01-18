clear all
close all

% Load coordinates of electrodes (from Susanne)
cap_coordinates_load = load('../data/Koordinaten_EEG2.mat');

% plot loaded 2D coordinates
if false
    figure
    hold on
    plot(cap_coordinates_load.M(:,1),cap_coordinates_load.M(:,2),'b*')
    plot(cap_coordinates_load.P(:,1),cap_coordinates_load.P(:,2),'r*')
    hold off
end
    
% define 2D coordinates of electrodes
cap_coordinates_2D = cap_coordinates_load.P;

% compute 3D coordinates of electrodes (magic)
cap_coordinates_3D = zeros(size(cap_coordinates_2D,1),3);
R = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
delta_t = 0.01;
projection_center = [-0.3,0,-3.0];
for k=1:size(cap_coordinates_3D,1)
    % scale
%    my_scale = [0.38 0.46];
    my_scale = [1 1];
    
    cap_coordinates_2D(k,:) = my_scale.*cap_coordinates_2D(k,:)*R + [-0.2 0];
    
    % x,y
    cap_coordinates_3D(k,:) = [cap_coordinates_2D(k,:), 2.3];
    direction = projection_center - cap_coordinates_3D(k,:);

    % z
    while norm((cap_coordinates_3D(k,:)+[0.0,0.0,0.7])./[3.5,2.9,3.0]) > 1 && cap_coordinates_3D(k,3) > -2.0
        cap_coordinates_3D(k,:) = cap_coordinates_3D(k,:) + delta_t*direction;
    end
end
    
if true
    % load brain data
    load('../data/brain_coordinates_simple.mat','fv');
    
    figure
    hold on

    [x,y,z] = sphere(30);
    x = 3.5*x;
    y = 2.9*y;
    z = 3.0*z;
    surf(x+0.13,y,z-0.7);     
    alpha 0.2;
    
    % plot brain
    patch(fv,'FaceColor',       [0.6 0.6 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
    
    % plot 2D
    plot3(cap_coordinates_2D(:,1),cap_coordinates_2D(:,2),2.3*ones(size(cap_coordinates_2D(:,1))),'r*')
 
    % plot 3D
    plot3(cap_coordinates_3D(:,1),cap_coordinates_3D(:,2),cap_coordinates_3D(:,3),'b*')
    for k=1:size(cap_coordinates_2D,1)
        text(cap_coordinates_3D(k,1),cap_coordinates_3D(k,2),cap_coordinates_3D(k,3),num2str(k));

        plot3([cap_coordinates_2D(k,1) projection_center(1)],[cap_coordinates_2D(k,2) projection_center(2)],[2.3 projection_center(3)],'b-')    
        
    end
    
    hold off

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');
    % Fix the axes scaling, and set a nice view angle
    axis('image');
%    view([-135 35]);
    view([-90 90]);

end

if true
   save('../data/cap_coordinates_3D.mat','cap_coordinates_3D');
end
