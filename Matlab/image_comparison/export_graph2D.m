% create petsc file with 2D graph example
%
% Lukas Pospisil, USI, Lugano 2016
%

addpath(genpath(fullfile(pwd,'../common')))

n = 2000; % number of points

ranges = [1 floor(n/4) floor(2/3*n) n];
mu{1} = [0;0];
mu{2} = [3;0]; 
mu{3} = [1.5; 2.5];
sigma{1} = 0.3*[1,0.5;0.5,1];
sigma{2} = 0.3*[3,0.5;0.5,1.5];
sigma{3} = 0.3*[3,0;0,3];

coordinates = zeros(2,n);
for j = 1:length(ranges)-1
    for i = ranges(j):ranges(j+1)
        % compute coordinate
        coordinates(:,i) = mvnrnd(mu{j},sigma{j},1)';
    end
end

% plot graph
figure
hold on
plot(coordinates(1,:),coordinates(2,:),'b*')
axis equal
hold off

% export coordinate to petscfile
savebin( 'output/test_graph2D.bin', [coordinates(1,:) coordinates(2,:)]' );
