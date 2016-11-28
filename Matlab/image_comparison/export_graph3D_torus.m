% create petsc file with 3D graph example
% points on torus
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all
addpath(genpath(fullfile(pwd,'../common')))

% number of points
a=5;
c=10;

[u,v]=meshgrid(0:5:360);
x=(c+a*cosd(v)).*cosd(u);
y=(c+a*cosd(v)).*sind(u);
z=a*sind(v);

coordinates = [];
for i=1:size(x,1)
    for j=1:size(x,2)
        coordinates(:,end+1) = [x(i,j);y(i,j);z(i,j)];
    end
end

% plot graph
figure
hold on
plot3(coordinates(1,:),coordinates(2,:),coordinates(3,:),'b*')
%surfl(x,y,z)
axis equal
hold off

% export coordinate to petscfile
savebin( 'output/test_graph3D_torus.bin', [coordinates(1,:) coordinates(2,:) coordinates(3,:)]' );
