function evaluationPointsInSkewCoordinates = computeSkewCoordinates(evaluationPointsInLocalCoordinates, nodalPoints, J0, shapeFunctionObject)
%COMPUTESKEWCOORDINATES Computes the skew coordinates for given local coordinates
%
%   Skew coordinates are used, e.g., in the context of Petrov-Galerkin formulations
%
%   REFERENCES
%   https://doi.org/10.1002/nme.6817
%   https://doi.org/10.1002/nme.7166
%
%   AUTHOR
%   Felix Zaehringer

% check & normalize input
assert(isa(shapeFunctionObject, 'shapeFunctionClass'), 'shapeFunctionObject must be of type shapeFunctionClass!');
dimension = shapeFunctionObject.dimension;
numberOfNodes = shapeFunctionObject.numberOfNodes;
assert(size(nodalPoints, 1) == dimension || size(nodalPoints, 2) == dimension, 'Size of nodalPoints matrix does not fit the given dimension!');
assert(size(evaluationPointsInLocalCoordinates, 1) == dimension || size(evaluationPointsInLocalCoordinates, 2) == dimension, 'Size of evaluationPointsInLocalCoordinates matrix does not fit the given dimension!');
if size(nodalPoints, 2) == dimension
    nodalPoints = nodalPoints.';
end
if size(evaluationPointsInLocalCoordinates, 2) == dimension
    evaluationPointsInLocalCoordinates = evaluationPointsInLocalCoordinates.';
end
assert(size(nodalPoints, 2) == numberOfNodes, 'Size of nodalPoints matrix does not fit the given numberOfNodes!');
assert(size(J0, 1) == dimension && size(J0, 2) == dimension, 'Size of J0 matrix does not fit the given dimension!')

% compute skew coordinates
numberOfEvaluationPoints = size(evaluationPointsInLocalCoordinates, 2);
N_k_I = computeLagrangeShapeFunction(dimension, numberOfNodes, numberOfEvaluationPoints, evaluationPointsInLocalCoordinates);
N0_k_I = shapeFunctionObject.N0_k_I;
evaluationPointsInSkewCoordinates = J0 \ (nodalPoints * (N_k_I - N0_k_I).');
end