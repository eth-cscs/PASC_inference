function [ brain_values ] = map_electrode_to_brain( brain_coordinates, cap_coordinates, cap_values )

% allocate vectors
brain_values = zeros(size(brain_coordinates,1),1); % output vector
alphas = zeros(size(cap_values,1),1); % coeficient of combination

% for each brain vertex
for i=1:length(brain_values)
    % coeficients are given by distance from electrodes
    for k=1:length(alphas)
       alphas(k) = 1/norm(brain_coordinates(i,:) - cap_coordinates(k,:));
    end
    
    % normalize coeficients
    alphas = alphas./sum(alphas);
    
    % compute new vertex value as a linear combination of values of all
    % electrodes
    brain_values(i) = sum(alphas.*cap_values);
end

end
