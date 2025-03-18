function evaluationPointsInQACMIICoordinates = computeQACMIICoordinates(evaluationPointsInLocalCoordinates, nodalPoints)
%COMPUTEQACMIICOORDINATES Computes the QACM-II coordinates for given local coordinates
%
%   QACM-II coordinates are used, e.g., in the context of Petrov-Galerkin formulations
%
%   REFERENCES
%   https://doi.org/10.1002/nme.6817
%   https://doi.org/10.1002/nme.7166
%
%   AUTHOR
%   Felix Zaehringer

% check & normalize input
if size(nodalPoints, 2) == 2
    nodalPoints = nodalPoints.';
end
assert(size(nodalPoints, 2) == 4, 'QACM-II coordinates are only defined for 4-node quadrilaterals!')
if size(evaluationPointsInLocalCoordinates, 2) == 2
    evaluationPointsInLocalCoordinates = evaluationPointsInLocalCoordinates.';
end

% compute QACM-II coordinates
% nodal coordinates
node1 = nodalPoints(:, 1);
node2 = nodalPoints(:, 2);
node3 = nodalPoints(:, 3);
node4 = nodalPoints(:, 4);
% diagonals
d1 = node3 - node1;
d2 = node4 - node2;
d1Length = norm(d1);
d2Length = norm(d2);
% side lengths
l12 = norm(node2-node1);
l23 = norm(node2-node3);
l14 = norm(node1-node4);
% compute element area A
theta = acos(d1.'*d2/(d1Length * d2Length));
A = 1 / 2 * d1Length * d2Length * sin(theta);
% compute area A123
s123 = 1 / 2 * (d1Length + l12 + l23);
A123 = sqrt(s123*(s123 - d1Length)*(s123 - l12)*(s123 - l23));
% compute area A124
s124 = 1 / 2 * (d2Length + l12 + l14);
A124 = sqrt(s124*(s124 - d2Length)*(s124 - l12)*(s124 - l14));
% shape parameters
g1 = (A123 - A124) / A;
g2 = (A - A124 - A123) / A;
evaluationPointsInQACMIICoordinates = evaluationPointsInLocalCoordinates + [g2; g1] .* evaluationPointsInLocalCoordinates(1, :) .* evaluationPointsInLocalCoordinates(2, :);
end