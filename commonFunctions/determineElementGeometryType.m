function elementGeometryType = determineElementGeometryType(dimension, numberOfNodes)
%DETERMINEELEMENTGEOMETRYTYPE Element geometry type.
%   This function determines the geometry type of an element based on it's
%   number of nodes and the dimension.
%
%   CALL
%   elementGeometryType = determineElementGeometryType(numberOfNodes, dimension)
%   numberOfNodes: number of nodes
%   dimension: dimension
%   elementGeometryType: the geometry type of the element
%
%   CREATOR(S)
%   Felix Zaehringer

if dimension == 1
    elementGeometryType = 'oneDimensional';
elseif dimension == 2
    if ismember(numberOfNodes, [3, 6])
        elementGeometryType = 'triangular';
    elseif ismember(numberOfNodes, [4, 8, 9])
        elementGeometryType = 'quadrilateral';
    else
        error('Element geometry type is not implemented!');
    end
elseif dimension == 3
    if ismember(numberOfNodes, [4, 10])
        elementGeometryType = 'tetrahedral';
    elseif ismember(numberOfNodes, [8, 20, 27])
        elementGeometryType = 'hexahedral';
    else
        error('Element geometry type is not implemented!');
    end
else
    error('Dimension not implemented!');
end
end
