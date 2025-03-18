function [nodes, edof, edofNeumann] = meshOneDimensional(lengthX, numberOfElementsX, order)
%MESHONEDIMENSIONAL Mesh for one dimensional domains
%   This function returns the nodes and the edof for one dimensional
%   domains.
%
%   CALL
%   [nodes, edof, edofNeumann] = meshOneDimensional(lengthX, numberOfElementsX, order)
%   lengthX: length of the body
%   numberOfElementsX: number of finite elements
%   order: order of the ansatz functions
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   edofNeumann: last node number. Can be used to apply Neumann boundary conditions
%
%   CREATOR(S)
%   Felix Zaehringer

%% check input
assert(lengthX > 0, 'lengthX must be positive!');
assert(mod(numberOfElementsX, 1) == 0 & numberOfElementsX > 0, 'numberOfElementsX must be an positive integer!');
assert(mod(order, 1) == 0, 'order must be an integer!')

%% compute nodes
nodes = linspace(-lengthX/2, lengthX/2, numberOfElementsX*order+1).';

%% compute edof
edof = zeros(numberOfElementsX, order+1);
for ii = 1:order + 1
    edof(:, ii) = ii:order:((numberOfElementsX - 1) * order) + ii;
end
edof = [edof(:, 1), edof(:, end), edof(:, 2:end-1)];

%% compute edofNeumann
edofNeumann = edof(end, 2);
end
