function [hourglassVector, hourglassFunction] = computeHourglassVectorAndFunction(dimension, localCoordinates)
%COMPUTEHOURGLASSVECTORANDFUNCTION Computes the hourglass vector and function
% 
%   REFERENCES
%   https://doi.org/10.1002/nme.6817
%
%   AUTHOR
%   Felix Zaehringer

xi = localCoordinates(1, :);
if dimension >= 2
    eta = localCoordinates(2, :);
end
if dimension >= 3
    zeta = localCoordinates(3, :);
end

if dimension == 2
    hourglassVector = [+1; -1; +1; -1];
    hourglassFunction = xi .* eta;
elseif dimension == 3
    hourglassVector = zeros(8, 4);
    hourglassVector(:, 1) = [+1; +1; -1; -1; -1; -1; +1; +1];
    hourglassVector(:, 2) = [+1; -1; -1; +1; -1; +1; +1; -1];
    hourglassVector(:, 3) = [+1; -1; +1; -1; +1; -1; +1; -1];
    hourglassVector(:, 4) = [-1; +1; -1; +1; +1; -1; +1; -1];

    hourglassFunction = zeros(size(localCoordinates, 2), 4);
    hourglassFunction(:, 1) = eta .* zeta;
    hourglassFunction(:, 2) = xi .* zeta;
    hourglassFunction(:, 3) = xi .* eta;
    hourglassFunction(:, 4) = xi .* eta .* zeta;
else
    error('Dimension not implemented!');
end
end