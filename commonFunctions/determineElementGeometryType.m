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
    if mod(sqrt(numberOfNodes), 1) == 0 || numberOfNodes == 8
        elementGeometryType = 'quadrilateral';
    elseif numberOfNodes == 2 %% for strings
        elementGeometryType = 'oneDimensional';
    elseif ismember(numberOfNodes, [3, 6])
        elementGeometryType = 'triangular';
    else
        error('Element geometry type is not implemented!');
    end
elseif dimension == 3
    if mod(nthroot(numberOfNodes, 3), 1) == 0 || numberOfNodes == 20
        elementGeometryType = 'hexahedral';
    elseif ismember(numberOfNodes, [4, 10])
        elementGeometryType = 'tetrahedral';
    elseif numberOfNodes == 2 %% for strings
        elementGeometryType = 'oneDimensional';
    else
        error('Element geometry type is not implemented!');
    end
else
    error('Dimension not implemented!');
end
end
