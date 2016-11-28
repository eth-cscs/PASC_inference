% create petsc file with 1D graph example
%
% Lukas Pospisil, USI, Lugano 2016
%

addpath(genpath(fullfile(pwd,'../common')))
n = 51; % size of graph

% compute coordinate
x_coordinate = zeros(n,1);
for i=1:floor(n/2)
    x_coordinate(floor(n/2)+1+i) = 0.01*i^2;
    x_coordinate(floor(n/2)+1-i) = 0.01*i^2;
end

% export coordinate to petscfile
savebin( 'output/test_graph1D.bin', x_coordinate );
