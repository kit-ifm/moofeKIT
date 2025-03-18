function [nodes, edof, bounEdof] = meshCurvedBeam(innerRadius, outerRadius, numberOfElements, order, serendipity)
%MESHCURVEDBEAM Mesh for the curved beam example
%   This function returns the nodes and the edof for the curved beam example
%
%   CALL
%   [nodes, edof, bounEdof] = meshTensionRodWithHole(radius, lengthX, lengthY, numberOfElementsFactor, order, serendipity)
%   innerRadius: inner radius of the curved beam
%   outerRadius: outer radius of the curved beam (thickness of the beam = outerRadius - innerRadius)
%   numberOfElements: number of elements
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions (true / false)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   bounEdof: Neumann boundary (top of the tension rod)
%
%   REFERENCE
%   -
%
%   CREATOR(S)
%   Felix Zaehringer

% check input
assert(innerRadius < outerRadius, 'innerRadius must be smaller than outerRadius!');

% compute nodal positions
[nodesOriginal, edof, bounEdof] = meshRectangle(1, 1, numberOfElements, 1, order, serendipity);
uniqueXValsOfNodes = unique(nodesOriginal(:, 1));
numberOfColumnsWithNodes = length(uniqueXValsOfNodes);
nodes = nodesOriginal;
angleValues = linspace(90, 0, numberOfColumnsWithNodes);
for ii = 1:numberOfColumnsWithNodes
    angle = angleValues(ii);
    xValOriginal = uniqueXValsOfNodes(ii);
    nodesinColumn = find(nodesOriginal(:, 1) == xValOriginal);
    startX = innerRadius * cosd(angle);
    endX = outerRadius * cosd(angle);
    startY = innerRadius * sind(angle);
    endY = outerRadius * sind(angle);
    nodes(nodesinColumn, 1) = linspace(startX, endX, length(nodesinColumn));
    nodes(nodesinColumn, 2) = linspace(startY, endY, length(nodesinColumn));
end

end
