clear all
close all

% Load brain mesh
load('../data/brain_coordinates.mat','fv');
load('../data/cap_coordinates_3D.mat','cap_coordinates_3D');

% add some data to vertices
fv.FaceVertexCData = zeros(size(fv.vertices,1),1);
for i = 1:length(fv.vertices)
   fv.FaceVertexCData(i) = fv.vertices(i,3); 
end

figure
hold on

% Render brain
patch('Faces', fv.faces,...
      'Vertices', fv.vertices, ...
      'FaceVertexCData', fv.FaceVertexCData,...
      'EdgeColor',       'none',        ...
      'FaceLighting',    'gouraud',     ...
      'AmbientStrength', 0.15, ...
      'FaceColor', 'interp', ...
      'CDataMapping', 'scaled');

% plot electrodes
plot3(cap_coordinates_3D(:,1),cap_coordinates_3D(:,2),cap_coordinates_3D(:,3),'k*');
for k=1:size(cap_coordinates_3D,1)
    text(cap_coordinates_3D(k,1),cap_coordinates_3D(k,2),cap_coordinates_3D(k,3),num2str(k));
end

% Add a camera light, and tone down the specular highlighting
%light;
%lighting phong;
%shading interp;
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');

colorbar

axis off
view([-135 35]);

