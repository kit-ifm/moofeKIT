function [nodes, edof, bounEdof] = meshThickWalledCylinder(order, serendipity)
%MESHCURVEDBEAM Mesh for the thick walled cylinder example
%   This function returns the nodes and the edof for the thick walled cylinder example
%   This example is commonly used to check elements for volumetric locking.
%
%   CALL
%   [nodes, edof, bounEdof] = meshThickWalledCylinder(order, serendipity)
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions (true / false)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   bounEdof: Neumann boundary (on the inside of the cylinder)
%
%   REFERENCE
%   https://doi.org/10.1016/0168-874X(85)90003-4
%
%   CREATOR(S)
%   Felix Zaehringer

[nodesOriginal, edof, bounEdofOriginal] = meshRectangle(1, 1, 9, 5, order, serendipity);

% compute nodal positions
nodesOriginalMainYValues = linspace(-0.5, 0.5, 6);
nodesMainRadialValues = [3, 3.5, 4.2, 5.2, 6.75, 9];
numberOfEdgeNodesPerElement = (length(find(nodesOriginal(:, 1) == 0.5)) - length(nodesOriginalMainYValues)) / 5;
nodes = nodesOriginal;
for ii = 1:length(nodesOriginalMainYValues)
    nodesOnLine = find(abs(nodesOriginal(:, 2)-nodesOriginalMainYValues(ii)) < 1e-8);
    angleValues = linspace(90, 0, length(nodesOnLine));
    xValues = nodesMainRadialValues(ii) * cosd(angleValues);
    yValues = nodesMainRadialValues(ii) * sind(angleValues);
    nodes(nodesOnLine, :) = [xValues.', yValues.'];
    if ii ~= 1 && numberOfEdgeNodesPerElement > 0 % compute nodal positions of edge nodes
        nodesOriginalEdgeYValues = linspace(nodesOriginalMainYValues(ii-1), nodesOriginalMainYValues(ii), numberOfEdgeNodesPerElement+2);
        nodesOriginalEdgeYValues = nodesOriginalEdgeYValues(2:end-1);
        nodesEdgeRadialValues = linspace(nodesMainRadialValues(ii-1), nodesMainRadialValues(ii), numberOfEdgeNodesPerElement+2);
        nodesEdgeRadialValues = nodesEdgeRadialValues(2:end-1);
        for jj = 1:numberOfEdgeNodesPerElement
            nodesOnLine = find(abs(nodesOriginal(:, 2)-nodesOriginalEdgeYValues(jj)) < 1e-8);
            angleValues = linspace(90, 0, length(nodesOnLine));
            xValues = nodesEdgeRadialValues(jj) * cosd(angleValues);
            yValues = nodesEdgeRadialValues(jj) * sind(angleValues);
            nodes(nodesOnLine, :) = [xValues.', yValues.'];
        end
    end
end

% bounEdof
bounEdof = bounEdofOriginal.SY1;

end
