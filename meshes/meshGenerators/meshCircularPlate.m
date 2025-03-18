function [nodes, edof, bounEdof] = meshCircularPlate(radius, numberOfElements, order, serendipity)
%MESHCIRCULARPLATE Mesh for the circular plate example
%   This function returns the nodes and the edof for the circular plate
%   example
%
%   CALL
%   [nodes, edof, bounEdof] = meshCircularPlate(radius, numberOfElements, order, serendipity)
%   radius: radius of the plate
%   numberOfElements: number of finite elements
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   bounEdof: Neumann boundary. Not implemented yet!
%
%   REFERENCE
%   https://doi.org/10.1115/1.3157679
%
%   CREATOR(S)
%   Felix Zaehringer

assert(order == 1, 'Not implemented yet!');
assert(serendipity == false, 'Not implemented yet!');
assert(numberOfElements > 0 && mod(sqrt(numberOfElements/3), 1) == 0, 'Invalid numberOfElements given!')

numberOfElementsPerPatch = numberOfElements / 3;
numberOfNodesPerPatchInOneDirection = sqrt(numberOfElementsPerPatch) * order + 1;
numberOfNodesPerPatch = (sqrt(numberOfElementsPerPatch) * order + 1)^2;

% compute nodal positions of center element patch
nodesCenterElementPatch = zeros(numberOfNodesPerPatch, 2);
for ii = 1:numberOfNodesPerPatchInOneDirection
    maxXVal = radius / 2 - (radius / 2 - radius / 2 * cosd(45)) * (ii - 1) / (numberOfNodesPerPatchInOneDirection - 1);
    maxYVal = radius / 2 - (radius / 2 - radius / 2 * cosd(45)) * (ii - 1) / (numberOfNodesPerPatchInOneDirection - 1);
    nodesCenterElementPatch((ii - 1)*numberOfNodesPerPatchInOneDirection+1:(ii - 1)*numberOfNodesPerPatchInOneDirection+numberOfNodesPerPatchInOneDirection, 1) = 0:maxXVal / (numberOfNodesPerPatchInOneDirection - 1):maxXVal;
    nodesCenterElementPatch(ii:numberOfNodesPerPatchInOneDirection:end, 2) = 0:maxYVal / (numberOfNodesPerPatchInOneDirection - 1):maxYVal;
end

% compute nodal positions of right element patch
nodesRightElementPatch = zeros(numberOfNodesPerPatch, 2);
for ii = 1:numberOfNodesPerPatchInOneDirection
    firstNodeX = nodesCenterElementPatch(numberOfNodesPerPatchInOneDirection*ii, 1);
    firstNodeY = nodesCenterElementPatch(numberOfNodesPerPatchInOneDirection*ii, 2);
    lastNodeX = cosd(45/(numberOfNodesPerPatchInOneDirection - 1)*(ii - 1)) * radius;
    lastNodeY = sind(45/(numberOfNodesPerPatchInOneDirection - 1)*(ii - 1)) * radius;
    nodesRightElementPatch((ii - 1)*numberOfNodesPerPatchInOneDirection+1:ii*numberOfNodesPerPatchInOneDirection, 1) = linspace(firstNodeX, lastNodeX, numberOfNodesPerPatchInOneDirection);
    nodesRightElementPatch((ii - 1)*numberOfNodesPerPatchInOneDirection+1:ii*numberOfNodesPerPatchInOneDirection, 2) = linspace(firstNodeY, lastNodeY, numberOfNodesPerPatchInOneDirection);
end

% compute nodal positions of top element patch
nodesTopElementPatch = zeros(numberOfNodesPerPatch, 2);
for ii = 1:numberOfNodesPerPatchInOneDirection
    firstNodeX = nodesCenterElementPatch(end-numberOfNodesPerPatchInOneDirection+ii, 1);
    firstNodeY = nodesCenterElementPatch(end-numberOfNodesPerPatchInOneDirection+ii, 2);
    lastNodeX = sind(45/(numberOfNodesPerPatchInOneDirection - 1)*(ii - 1)) * radius;
    lastNodeY = cosd(45/(numberOfNodesPerPatchInOneDirection - 1)*(ii - 1)) * radius;
    nodesTopElementPatch(ii:numberOfNodesPerPatchInOneDirection:numberOfNodesPerPatchInOneDirection*(numberOfNodesPerPatchInOneDirection - 1)+ii, 1) = linspace(firstNodeX, lastNodeX, numberOfNodesPerPatchInOneDirection);
    nodesTopElementPatch(ii:numberOfNodesPerPatchInOneDirection:numberOfNodesPerPatchInOneDirection*(numberOfNodesPerPatchInOneDirection - 1)+ii, 2) = linspace(firstNodeY, lastNodeY, numberOfNodesPerPatchInOneDirection);
end

% create edof for element patches
[~, edofCenterElementPatch, ~] = meshRectangle(radius/2, radius/2, sqrt(numberOfElementsPerPatch), sqrt(numberOfElementsPerPatch), order, serendipity);
[~, edofRightElementPatch, ~] = meshRectangle(radius/2, radius/2, sqrt(numberOfElementsPerPatch), sqrt(numberOfElementsPerPatch), order, serendipity);
[~, edofTopElementPatch, ~] = meshRectangle(radius/2, radius/2, sqrt(numberOfElementsPerPatch), sqrt(numberOfElementsPerPatch), order, serendipity);

% tie mesh
edofRightElementPatch = edofRightElementPatch + max(max(edofCenterElementPatch));
[nodes, edof] = meshTie(nodesCenterElementPatch, edofCenterElementPatch, nodesRightElementPatch, edofRightElementPatch, eps);
edofTopElementPatch = edofTopElementPatch + max(max(edof));
[nodes, edof] = meshTie(nodes, edof, nodesTopElementPatch, edofTopElementPatch, eps);

bounEdof = struct();
end
