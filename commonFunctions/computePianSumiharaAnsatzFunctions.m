function [M, dM_xi] = computePianSumiharaAnsatzFunctions(continuumObject, dimension, orderAnsatzFunction, numberOfGausspoints, gaussPoints)
%COMPUTEPIANSUMIHARAANSATZFUNCTION Computes the ansatz functions for Pian-Sumihara elements.
%
%   CALL
%   [M, dM_xi] = computePianSumiharaAnsatzFunctions(continuumObject, dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
%   continuumObject: object of solidSuperClass
%   dimension: dimension
%   numberOfNodes: number of nodes
%   numberOfGausspoints: number of Gauss points
%   gaussPoint: Gauss points
%   M: Ansatz function
%   dM_xi:  Derivatives of the ansatz function with respect to the
%           local coordinates
%
%   CREATOR(S)
%   Felix Zaehringer

M = zeros((dimension^2 + dimension)/2*numberOfGausspoints, orderAnsatzFunction);
dM_xi = [];

stepSize = (dimension^2 + dimension) / 2;

xi = gaussPoints(1, :);
eta = gaussPoints(2, :);

if dimension == 2
    if orderAnsatzFunction == 5
        M(1:stepSize:end, 1) = 1;
        M(2:stepSize:end, 2) = 1;
        M(3:stepSize:end, 3) = 1;
        M(1:stepSize:end, 4) = eta;
        M(2:stepSize:end, 5) = xi;
    elseif orderAnsatzFunction == 8
        M(1:stepSize:end, 1) = 1;
        M(2:stepSize:end, 2) = 1;
        M(3:stepSize:end, 3) = 1;
        M(1:stepSize:end, 4) = eta;
        M(2:stepSize:end, 5) = xi;
        M(1:stepSize:end, 6) = xi .* eta;
        M(2:stepSize:end, 7) = xi .* eta;
        M(3:stepSize:end, 8) = xi .* eta;
    elseif orderAnsatzFunction == 12
        M(1:stepSize:end, 1) = 1;
        M(2:stepSize:end, 2) = 1;
        M(3:stepSize:end, 3) = 1;
        M(1:stepSize:end, 4) = xi;
        M(1:stepSize:end, 5) = eta;
        M(1:stepSize:end, 6) = xi .* eta;
        M(2:stepSize:end, 7) = xi;
        M(2:stepSize:end, 8) = eta;
        M(2:stepSize:end, 9) = xi .* eta;
        M(3:stepSize:end, 10) = xi;
        M(3:stepSize:end, 11) = eta;
        M(3:stepSize:end, 12) = xi .* eta;
    else
        error('Ansatz function not implemented!');
    end
elseif dimension == 3
    zeta = gaussPoints(3, :);

    if orderAnsatzFunction == 18
        M(1:stepSize:end, 1) = 1;
        M(2:stepSize:end, 2) = 1;
        M(3:stepSize:end, 3) = 1;
        M(4:stepSize:end, 4) = 1;
        M(5:stepSize:end, 5) = 1;
        M(6:stepSize:end, 6) = 1;
        M(1:stepSize:end, 7) = eta;
        M(1:stepSize:end, 8) = zeta;
        M(2:stepSize:end, 9) = zeta;
        M(2:stepSize:end, 10) = xi;
        M(3:stepSize:end, 11) = xi;
        M(3:stepSize:end, 12) = eta;
        M(4:stepSize:end, 13) = zeta;
        M(5:stepSize:end, 14) = xi;
        M(6:stepSize:end, 15) = eta;
        M(1:stepSize:end, 16) = eta .* zeta;
        M(2:stepSize:end, 17) = zeta .* xi;
        M(3:stepSize:end, 18) = xi .* eta;
    else
        error('Ansatz function not implemented!');
    end
else
    error('No ansatz functions implemented for this dimension!');
end

end
