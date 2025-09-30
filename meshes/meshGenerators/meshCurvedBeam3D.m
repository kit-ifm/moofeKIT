function [nodes, edof, bounEdof] = meshCurvedBeam3D(innerRadius, outerRadius, numberOfElements,nelZ, width, order, serendipity)
%MESHCURVEDBEAM3D Mesh for the curved beam example
%   This function returns the nodes and the edof for the curved beam example
%
%   CALL
%   [nodes, edof, edofneumann] = meshCurvedBeam3D(innerRadius, outerRadius, numberOfElements, nelZ, width, order, serendipity)
%   innerRadius: inner radius of the curved beam
%   outerRadius: outer radius of the curved beam (thickness of the beam = outerRadius - innerRadius)
%   numberOfElements: number of elements in x/y-plane
%   nelZ: number of elements in z-Dimension
%   width: width of the beam
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions (true / false)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   edofNeumann: Neumann boundary 
%   
%   REFERENCE
%   -
%
%   SEE ALSO /based on
%   meshCurvedBeam 
%   meshGeneratorCube
%
%   CREATOR(S)
%   Jakob Hammes

% check input
assert(innerRadius < outerRadius, 'innerRadius must be smaller than outerRadius!');
assert(numberOfElements >= 1, 'number of elements must be positve and greater or equal to 1!');
assert(nelZ >= 1, 'number of elements in z-Dimension must be greater or equal to 1');

% compute nodal positions 
[nodesOriginal, edof, bounEdof] = meshGeneratorCube(1, 1,width, numberOfElements, 1,nelZ, order, serendipity);

uniqueXValuesOfNodes = unique(nodesOriginal(:, 1));
angleValues = linspace(90, 0, (numberOfElements+1));
radii = [innerRadius outerRadius];
nodes = nodesOriginal;

for ii = 1:(numberOfElements+1)
    angle = angleValues(ii);
    xValueOriginal = uniqueXValuesOfNodes(ii);
    nodesInColumn = find(nodesOriginal(:,1)== xValueOriginal);
    xValues = radii.*cosd(angle);
    yValues = radii.*sind(angle);
    

    nodes(nodesInColumn(1:2),1) = xValues;
    nodes(nodesInColumn(1:2),2) = yValues;
end

numberOfnodesInPlane=size(nodes,1)/(nelZ+1);
for ii = 1 : nelZ
    nodes(numberOfnodesInPlane*ii+1:numberOfnodesInPlane*(ii+1),1:2) = nodes(1:numberOfnodesInPlane,1:2);
end
