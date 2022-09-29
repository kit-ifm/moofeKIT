function [N, dNr] = computePianSumiharaAnsatzFunctions(dimension, orderAnsatzFunction, numberOfGausspoints, gaussPoints)
%COMPUTEPIANSUMIHARAANSATZFUNCTION Computes the ansatz functions for Pian-Sumihara elements.
%
%   CALL
%   [N, dNr] = computePianSumiharaAnsatzFunctions(dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
%   dimension: dimension
%   numberOfNodes: number of nodes
%   numberOfGausspoints: number of Gauss points
%   gaussPoint: Gauss points
%   N: Ansatz function
%   dNr: Derivatives of the ansatz function with respect to the
%        isoparametric coordinates
%
%   CREATOR(S)
%   Felix Zaehringer

N = zeros((dimension^2 + dimension)/2*numberOfGausspoints, orderAnsatzFunction);
dNr = [];

stepSize = (dimension^2 + dimension) / 2;

if dimension == 2
    if orderAnsatzFunction == 5
        N(1:stepSize:end, 1) = 1;
        N(2:stepSize:end, 2) = 1;
        N(3:stepSize:end, 3) = 1;
        N(1:stepSize:end, 4) = gaussPoints(:, 2);
        N(2:stepSize:end, 5) = gaussPoints(:, 1);
    elseif orderAnsatzFunction == 8
        N(1:stepSize:end, 1) = 1;
        N(2:stepSize:end, 2) = 1;
        N(3:stepSize:end, 3) = 1;
        N(1:stepSize:end, 4) = gaussPoints(:, 2);
        N(2:stepSize:end, 5) = gaussPoints(:, 1);
        N(1:stepSize:end, 6) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(2:stepSize:end, 7) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(3:stepSize:end, 8) = gaussPoints(:, 1) .* gaussPoints(:, 2);
    elseif orderAnsatzFunction == 12
        N(1:stepSize:end, 1) = 1;
        N(2:stepSize:end, 2) = 1;
        N(3:stepSize:end, 3) = 1;
        N(1:stepSize:end, 4) = gaussPoints(:, 1);
        N(1:stepSize:end, 5) = gaussPoints(:, 2);
        N(1:stepSize:end, 6) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(2:stepSize:end, 7) = gaussPoints(:, 1);
        N(2:stepSize:end, 8) = gaussPoints(:, 2);
        N(2:stepSize:end, 9) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(3:stepSize:end, 10) = gaussPoints(:, 1);
        N(3:stepSize:end, 11) = gaussPoints(:, 2);
        N(3:stepSize:end, 12) = gaussPoints(:, 1) .* gaussPoints(:, 2);
    else
        error('Ansatz function not implemented!');
    end
else
    error('No ansatz functions implemented for this dimension!');
end

end
