function [nodes] = adaptMeshToCornerNodes(nodes, edof, dimension)
%ADAPTMESHTOCORNERNODES Adjust the entire mesh to (distorted) corner nodes
%   This function adjusts the nodal coordinates of the entire FE mesh, such
%   that it fits to the coordinates of the corner nodes. More precisely,
%   this means that you can distort the corner nodes of an mesh and adjust
%   the edge + inner nodes with this function.
%
%   INPUT:
%   nodes: nodal coordinates
%   edof: edof
%
%   OUTPUT:
%   nodes
%
%   CREATOR(S):
%   Felix Zaehringer

numberOfElements = size(edof, 1);
numberOfNodesPerElement = size(edof, 2);

% check input
assert(dimension == 2, 'This function only works for 2D yet!');
elementGeometryType = determineElementGeometryType(dimension, numberOfNodesPerElement);
assert(strcmpi(elementGeometryType, 'quadrilateral'), 'This function only works for quadrilaterals yet!');
assert(numberOfNodesPerElement <= 9, 'This function only works for quadrilaterals with up to 9 nodes yet!')

% determine additional information
isSerendipity = true;
if mod(sqrt(numberOfNodesPerElement), 1) == 0
    isSerendipity = false;
end

% adjust node positions
for e=1:numberOfElements
    localEdof = edof(e, :);
    cornerNodes = nodes(localEdof(1:4), :);
    if numberOfNodesPerElement > 4
        nodes(localEdof(5), :) = (cornerNodes(1, :) + cornerNodes(2, :)) / 2;
        nodes(localEdof(6), :) = (cornerNodes(2, :) + cornerNodes(3, :)) / 2;
        nodes(localEdof(7), :) = (cornerNodes(3, :) + cornerNodes(4, :)) / 2;
        nodes(localEdof(8), :) = (cornerNodes(1, :) + cornerNodes(4, :)) / 2;
        if ~isSerendipity
            nodes(localEdof(9), :) = (cornerNodes(1, :) + cornerNodes(2, :) + cornerNodes(3, :) + cornerNodes(4, :)) / 4;
        end
    end
end
