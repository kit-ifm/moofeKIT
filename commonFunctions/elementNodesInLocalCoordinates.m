function nodesXi = elementNodesInLocalCoordinates(dimension, elementGeometryType, elementNumberOfNodes)
% ELEMENTNODESINLOCALCOORDINATES Local coordinates of the reference element
% Returns the nodes of the reference element in xi-coordinate
% frame. Node numbering has to fit lagrange shape functions of moofeKIT.
%
% CREATOR(S)
% Robin Pfefferkorn, Felix Zaehringer

nodesXi = zeros(dimension, elementNumberOfNodes);

if dimension == 1
    if strcmpi(elementGeometryType, 'oneDimensional')
        if elementNumberOfNodes == 1
            nodesXi = 0;
        elseif elementNumberOfNodes > 1 && mod(elementNumberOfNodes, 1) == 0
            [nodes, edof] = meshOneDimensional(2, 1, elementNumberOfNodes-1);
            nodesXi = nodes(edof)';
        else
            error('Element number of nodes not implemented');
        end
    else
        error('Element geometry type not implemented');
    end
elseif dimension == 2
    if strcmpi(elementGeometryType, 'triangular')
        if elementNumberOfNodes == 3
            % 3-nodes
            nodesXi(1, :) = [+0, +1, +0];
            nodesXi(2, :) = [+0, +0, +1];
        elseif elementNumberOfNodes == 6
            % 6-nodes
            nodesXi(1, :) = [+0, +1, +0, +1 / 2, +1 / 2, +0];
            nodesXi(2, :) = [+0, +0, +1, +0, +1 / 2, +1 / 2];
        else
            error('Element number of nodes not implemented');
        end
    elseif strcmpi(elementGeometryType, 'quadrilateral')
        if mod(sqrt(elementNumberOfNodes), 1) == 0
            % lagrange elements
            [nodes, edof] = meshRectangle(2, 2, 1, 1, sqrt(elementNumberOfNodes)-1);
            nodesXi = nodes(edof, :)';
        elseif elementNumberOfNodes == 8
            % 8-node serendipity element
            [nodes, edof] = meshRectangle(2, 2, 1, 1, 2);
            nodesXi = nodes(edof, :)';
            nodesXi = nodesXi(:, 1:end-1);
        else
            error('Element number of nodes not implemented');
        end
    else
        error('Element geometry type not implemented');
    end
elseif dimension == 3
    if strcmpi(elementGeometryType, 'tetrahedral')
        if elementNumberOfNodes == 4
            % 4-node
            nodesXi(1, :) = [+0, +1, +0, +0];
            nodesXi(2, :) = [+0, +0, +1, +0];
            nodesXi(3, :) = [+0, +0, +0, +1];
        elseif elementNumberOfNodes == 10
            % 10-nodes
            nodesXi(1, :) = [+0, +1, +0, +0, 1 / 2, 1 / 2, 0, 0, 1 / 2, 0];
            nodesXi(2, :) = [+0, +0, +1, +0, 0, 1 / 2, 1 / 2, 0, 0, 1 / 2];
            nodesXi(3, :) = [+0, +0, +0, +1, 0, 0, 0, 1 / 2, 1 / 2, 1 / 2];
        else
            error('Element number of nodes not implemented');
        end
    elseif strcmpi(elementGeometryType, 'hexahedral')
        if elementNumberOfNodes == 8
            % 8-node tri linear lagrange element
            [nodes, edof] = meshGeneratorCube(2, 2, 2, 1, 1, 1, 1, false);
            nodesXi = nodes(edof, :)';
        elseif elementNumberOfNodes == 20
            % 20-node serendipity element
            [nodes, edof] = meshGeneratorCube(2, 2, 2, 1, 1, 1, 2, true);
            nodesXi = nodes(edof, :)';
        elseif elementNumberOfNodes == 27
            % 27-node tri quadratic lagrange element
            [nodes, edof] = meshGeneratorCube(2, 2, 2, 1, 1, 1, 2, false);
            nodesXi = nodes(edof, :)';
        else
            error('Element number of nodes not implemented');
        end
    else
        error('Element geometry type not implemented');
    end
else
    error('Element nodes not implemented for this dimension!');
end
end