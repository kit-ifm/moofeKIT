function [N, dNr] = computeEASAnsatzFunctions(dimension, orderAnsatzFunction, numberOfGausspoints, gaussPoints)
%COMPUTEEASANSATZFUNCTION Computes the ansatz functions for EAS elements.
%
%   CALL
%   [N, dNr] = computeEASAnsatzFunctions(dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
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
    if orderAnsatzFunction == 4
        N(1:stepSize:end, 1) = gaussPoints(:, 1);
        N(2:stepSize:end, 2) = gaussPoints(:, 2);
        N(3:stepSize:end, 3) = gaussPoints(:, 1);
        N(3:stepSize:end, 4) = gaussPoints(:, 2);
    elseif orderAnsatzFunction == 7
        N(1:stepSize:end, 1) = gaussPoints(:, 1);
        N(2:stepSize:end, 2) = gaussPoints(:, 2);
        N(3:stepSize:end, 3) = gaussPoints(:, 1);
        N(3:stepSize:end, 4) = gaussPoints(:, 2);
        N(1:stepSize:end, 5) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(2:stepSize:end, 6) = gaussPoints(:, 1) .* gaussPoints(:, 2);
        N(3:stepSize:end, 7) = gaussPoints(:, 1) .* gaussPoints(:, 2);
    else
        error('Ansatz function not implemented!');
    end
elseif dimension == 3
    if orderAnsatzFunction == 9
        N(1:stepSize:end, 1) = gaussPoints(:, 1);
        N(2:stepSize:end, 2) = gaussPoints(:, 2);
        N(3:stepSize:end, 3) = gaussPoints(:, 3);
        N(4:stepSize:end, 4) = gaussPoints(:, 1);
        N(4:stepSize:end, 5) = gaussPoints(:, 2);
        N(5:stepSize:end, 6) = gaussPoints(:, 2);
        N(5:stepSize:end, 7) = gaussPoints(:, 3);
        N(6:stepSize:end, 8) = gaussPoints(:, 1);
        N(6:stepSize:end, 9) = gaussPoints(:, 3);
    elseif orderAnsatzFunction == 12
        N(1:stepSize:end, 1) = gaussPoints(:, 1);
        N(2:stepSize:end, 2) = gaussPoints(:, 2);
        N(3:stepSize:end, 3) = gaussPoints(:, 3);
        N(4:stepSize:end, 4) = gaussPoints(:, 1);
        N(4:stepSize:end, 5) = gaussPoints(:, 2);
        N(5:stepSize:end, 6) = gaussPoints(:, 2);
        N(5:stepSize:end, 7) = gaussPoints(:, 3);
        N(6:stepSize:end, 8) = gaussPoints(:, 1);
        N(6:stepSize:end, 9) = gaussPoints(:, 3);
        N(1:stepSize:end, 10:12) = [gaussPoints(:, 1) .* gaussPoints(:, 2), gaussPoints(:, 2) .* gaussPoints(:, 3), gaussPoints(:, 3) .* gaussPoints(:, 1)];
        N(2:stepSize:end, 10:12) = [gaussPoints(:, 1) .* gaussPoints(:, 2), gaussPoints(:, 2) .* gaussPoints(:, 3), gaussPoints(:, 3) .* gaussPoints(:, 1)];
        N(3:stepSize:end, 10:12) = [gaussPoints(:, 1) .* gaussPoints(:, 2), gaussPoints(:, 2) .* gaussPoints(:, 3), gaussPoints(:, 3) .* gaussPoints(:, 1)];
    else
        error('Ansatz function not implemented!');
    end
else
    error('No ansatz functions implemented for this dimension!');
end

end
