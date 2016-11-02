% create petsc file with 3D graph example
% some random points on sphere
% (maybe real cities on Earth?)
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all
addpath('common')
n = 3000; % number of points
radius = 1; % radius of the sphere

ranges = [1 floor(n/4) floor(3/8*n) floor(n/2) floor(3/5*n) floor(5/7*n) floor(5/6*n) n];
mu{1} = [-2;-2;-2];
mu{2} = [0;3;1.5]; 
mu{3} = [-1.5; 2.5; 0];
mu{4} = [-4;4;0]; 
mu{5} = [3;-4.5;0]; 
mu{6} = [2;0;0]; 
mu{7} = [0;0;5]; 

sigma{1} = 0.3*[1,0.5,0; 0.5, 1,0.5; 0 0.5, 5];
sigma{2} = 0.5*[5,0.5,0; 0.5, 1,0.5; 0 0.5, 1];
sigma{3} = 0.6*[3,0,0;0,3,0;0,0,3];
sigma{4} = 0.7*[3,0,0.5;0,5,0.2;0.5,0.2,3];
sigma{5} = 1.6*sigma{2};
sigma{6} = 1.0*sigma{1};
sigma{7} = 1.3*sigma{3};

coordinates = zeros(3,n);
for j = 1:length(ranges)-1
    for i = ranges(j):ranges(j+1)
        % compute coordinate
        coordinates(:,i) = mvnrnd(mu{j},sigma{j},1)';
        
        % project onto sphere
        r2 = norm(coordinates(:,i));
        coordinates(:,i) = radius*coordinates(:,i)/r2;
    end
end

% plot graph
figure
hold on
plot3(coordinates(1,:),coordinates(2,:),coordinates(3,:),'b*')
%surfl(x,y,z)
xlabel('x1')
ylabel('x2')
zlabel('x3')
axis equal
hold off

% export coordinate to petscfile
savebin( 'output/test_graph3D_sphere.bin', [coordinates(1,:) coordinates(2,:) coordinates(3,:)]' );
