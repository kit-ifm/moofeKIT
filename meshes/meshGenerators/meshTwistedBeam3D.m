function [nodes, edof, bounEdof] = meshTwistedBeam3D(length, height, width, angle, nelX, nelY, nelZ, order, serendipity)
%MESHTWISTEDBEAM3D Mesh for the twisted beam example
%   This function returns the nodes and the edof for the twisted beam example
%
%   CALL
%   [nodes, edof, edofneumann] = meshTwistedBeam3D(length, height, width, angle, nelX, nelY, nelZ, order, serendipity)
%   lenght: lenght of the beam
%   height: height of the beam
%   width:  width of the beam
%   angle:  angle of rotation at the end of the beam
%   nelX: number of elements in x-Dimension
%   nelY: number of elements in y-Dimension
%   nelZ: number of elements in z-Dimension
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions (true / false)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   bounEdof: Neumann boundary 
%   
%   REFERENCE
%   -
%
%   SEE ALSO /based on
%   meshGeneratorCube
%   meshCurvedBeam3D
%
%   CREATOR(S)
%   Jakob Hammes

% check input
assert(nelX >= 1, 'number of elements in x-Dimension must be greater or equal to 1');
assert(nelY >= 1, 'number of elements in y-Dimension must be greater or equal to 1');
assert(nelZ >= 1, 'number of elements in z-Dimension must be greater or equal to 1');

% compute nodal positions of undistorted beam
[nodesOriginal, edof, bounEdof] = meshGeneratorCube(length, height, width, nelX, nelY, nelZ, order, serendipity);

% shift nodes 
nodesOriginal(:,1) = nodesOriginal(:,1)+length/2;

% rotate cube around x-Axis
nodes = zeros(size(nodesOriginal));
uniqueXValuesOfNodes = unique(nodesOriginal(:, 1));
angleValues = linspace(0, angle, (nelX+1));

 for ii = 1:(nelX+1)
     currentAngle = angleValues(ii);
     xValueOriginal = uniqueXValuesOfNodes(ii);
     nodesInColumn = find(nodesOriginal(:,1)== xValueOriginal);
     nodes(nodesInColumn,1) = nodesOriginal(nodesInColumn,1);
     nodes(nodesInColumn,2:3) = ([cosd(currentAngle) -sind(currentAngle); sind(currentAngle) cosd(currentAngle)]*nodesOriginal(nodesInColumn,2:3)')';
 end


