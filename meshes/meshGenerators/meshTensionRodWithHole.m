function [nodes, edof, bounEdof] = meshTensionRodWithHole(radius, lengthX, lengthY, numberOfElementsFactor, order, serendipity)
%MESHTENSIONRODWITHHOLE Mesh for the tension rod with a hole example
%   This function returns the nodes and the edof for the tension rod with a
%   hole example
%
%   CALL
%   [nodes, edof, bounEdof] = meshTensionRodWithHole(radius, lengthX, lengthY, numberOfElementsFactor, order, serendipity)
%   radius: radius of the hole
%   lengthX: width of the tension rod (including the hole)
%   lengthY: height of the tension rod (including the hole)
%   numberOfElementsFactor: totalNumberOfElements = 3*(numberOfElementsFactor)^2
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
assert(radius < lengthX, 'radius must be smaller than lengthX!');
assert(lengthX < lengthY, 'lengthX must be smaller than lengthY!');

% compute nodal positions of bottom right element patch
[nodesBottomRightElementPatchOriginal, edofBottomRightElementPatch] = meshRectangle(1, 1, numberOfElementsFactor, numberOfElementsFactor, order, serendipity);
uniqueYValsOfNodes = unique(nodesBottomRightElementPatchOriginal(:, 2));
numberOfColumnsWithNodes = length(uniqueYValsOfNodes);
nodesBottomRightElementPatch = nodesBottomRightElementPatchOriginal;
for ii = 1:numberOfColumnsWithNodes
    angle = atand((ii - 1)/(numberOfColumnsWithNodes - 1));
    yValOriginal = uniqueYValsOfNodes(ii);
    nodesinRow = find(nodesBottomRightElementPatchOriginal(:, 2) == yValOriginal);
    startX = radius * cosd(angle);
    endX = lengthX;
    startY = radius * sind(angle);
    endY = lengthX * (ii - 1) / (numberOfColumnsWithNodes - 1);
    nodesBottomRightElementPatch(nodesinRow, 1) = linspace(startX, endX, length(nodesinRow));
    nodesBottomRightElementPatch(nodesinRow, 2) = linspace(startY, endY, length(nodesinRow));
end

% compute nodal positions of middle element patch
[nodesMiddleElementPatchOriginal, edofMiddleElementPatch] = meshRectangle(1, 1, numberOfElementsFactor, numberOfElementsFactor, order, serendipity);
uniqueXValsOfNodes = unique(nodesMiddleElementPatchOriginal(:, 1));
numberOfColumnsWithNodes = length(uniqueXValsOfNodes);
nodesMiddleElementPatch = nodesMiddleElementPatchOriginal;
for ii = 1:numberOfColumnsWithNodes
    angle = atand((ii - 1)/(numberOfColumnsWithNodes - 1));
    xValOriginal = uniqueXValsOfNodes(ii);
    nodesinRow = find(nodesMiddleElementPatchOriginal(:, 1) == xValOriginal);
    startX = radius * sind(angle);
    endX = lengthX * (ii - 1) / (numberOfColumnsWithNodes - 1);
    startY = radius * cosd(angle);
    endY = lengthX;
    nodesMiddleElementPatch(nodesinRow, 1) = linspace(startX, endX, length(nodesinRow));
    nodesMiddleElementPatch(nodesinRow, 2) = linspace(startY, endY, length(nodesinRow));
end

% compute nodal positions of top element patch
[nodesTopElementPatch, edofTopElementPatch, bounEdofTopElementPatch] = meshRectangle(lengthX, lengthY-lengthX, numberOfElementsFactor, numberOfElementsFactor, order, serendipity);
nodesTopElementPatch(:, 1) = nodesTopElementPatch(:, 1) + lengthX / 2;
nodesTopElementPatch(:, 2) = nodesTopElementPatch(:, 2) + (lengthY - lengthX) / 2 + lengthX;

% tie mesh
edofMiddleElementPatch = edofMiddleElementPatch + max(max(edofTopElementPatch));
[nodes, edof] = meshTie(nodesMiddleElementPatch, edofMiddleElementPatch, nodesTopElementPatch, edofTopElementPatch, 1e-8);
edofBottomRightElementPatch = edofBottomRightElementPatch + max(max(edof));
[nodes, edof] = meshTie(nodes, edof, nodesBottomRightElementPatch, edofBottomRightElementPatch, 1e-8);

% bounEdof
nodesForXVals = find(nodesMiddleElementPatch(:, 2) == lengthX);
numberOfBoundaryElements = numberOfElementsFactor;
bounEdof = zeros(size(bounEdofTopElementPatch.SY2));
for e = 1:numberOfBoundaryElements
    numberOfNodes = size(bounEdof, 2);
    index = 1 + (numberOfNodes - 1) * (e - 1):1 + (numberOfNodes - 1) * e;
    nodesForElement = nodesMiddleElementPatch(nodesForXVals(index), 1);
    for ii = 1:length(nodesForElement)
        bounEdof(e, ii) = find(abs(nodes(:, 2)-lengthY) < 1e-8 & abs(nodes(:, 1)-nodesForElement(ii)) < 1e-8);
    end
end
bounEdof = [bounEdof(:, 1), bounEdof(:, end), bounEdof(:, 2:end-1)];
end
