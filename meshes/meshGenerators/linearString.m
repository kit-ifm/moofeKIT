function [nodes, edof, edofNeumann] = linearString(length, numberOfElements, order, startPoint, endPoint)
%LINEARSTRING2D Mesh for one dimensional string domain
%   This function returns the nodes and the edof for one dimensional
%   string domains with two DOF at each node (i.e. 2D space).
%
%   CALL
%   [nodes, edof, edofNeumann] = linearMeshString2D(length, numberOfElements, order, endPoint)
%   length: length of the string
%   numberOfElements: number of finite elements
%   order: order of the ansatz functions
%   endPoint: Position in R2 of the right end, left end situated in origin (0,0)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   edofNeumann: last node number. Can be used to apply Neumann boundary conditions
%
%   CREATOR(S)
%   Felix Zaehringer, Philipp Kinon

%% check input
assert(length > 0, 'lengthX must be positive!');
assert(mod(numberOfElements, 1) == 0 & numberOfElements > 0, 'numberOfElementsX must be an positive integer!');
assert(mod(order, 1) == 0, 'order must be an integer!')

%% compute nodes
mat_nodes = (0 : length / (numberOfElements * order):length)';
dimension = max(size(endPoint));
nodes = zeros(numberOfElements * order +1, dimension);
for k = 1:dimension
    nodes(:,k) = startPoint(k) + mat_nodes*(endPoint(k)-startPoint(k))/length;
end
%% compute edof
edof = zeros(numberOfElements, order+1);
for ii = 1:order + 1
    edof(:, ii) = ii:order:((numberOfElements - 1) * order) + ii;
end
edof = [edof(:, 1), edof(:, end), edof(:, 2:end-1)];

%% compute edofNeumann
edofNeumann = edof(end, 2);
end
