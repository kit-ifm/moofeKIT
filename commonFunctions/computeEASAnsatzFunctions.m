function [M, dM_xi] = computeEASAnsatzFunctions(continuumObject, dimension, orderAnsatzFunction, numberOfGausspoints, gaussPoints)
%COMPUTEEASANSATZFUNCTION Computes the ansatz functions for EAS elements.
%
%   CALL
%   [M, dM_xi] = computeEASAnsatzFunctions(continuumObject, dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
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

xi = gaussPoints(1, :);
eta = gaussPoints(2, :);

if isa(continuumObject, 'plateClass') && contains(continuumObject.elementNameAdditionalSpecification, 'SimoRifai', IgnoreCase = true)
    M = zeros(2*numberOfGausspoints, orderAnsatzFunction);
    dM_xi = [];

    stepSize = 2;

    if orderAnsatzFunction == 4
        M(1:stepSize:end, 1) = xi;
        M(2:stepSize:end, 2) = eta;
        M(1:stepSize:end, 3) = xi .* eta;
        M(2:stepSize:end, 4) = xi .* eta;
    else
        error('Ansatz function not implemented!');
    end
else
    M = zeros((dimension^2 + dimension)/2*numberOfGausspoints, orderAnsatzFunction);
    dM_xi = [];

    stepSize = (dimension^2 + dimension) / 2;

    if dimension == 2
        if orderAnsatzFunction == 4
            M(1:stepSize:end, 1) = xi;
            M(2:stepSize:end, 2) = eta;
            M(3:stepSize:end, 3) = xi;
            M(3:stepSize:end, 4) = eta;
        elseif orderAnsatzFunction == 7
            M(1:stepSize:end, 1) = xi;
            M(2:stepSize:end, 2) = eta;
            M(3:stepSize:end, 3) = xi;
            M(3:stepSize:end, 4) = eta;
            M(1:stepSize:end, 5) = xi .* eta;
            M(2:stepSize:end, 6) = xi .* eta;
            M(3:stepSize:end, 7) = xi .* eta;
        else
            error('Ansatz function not implemented!');
        end
    elseif dimension == 3
        zeta = gaussPoints(3, :);

        if orderAnsatzFunction == 9
            M(1:stepSize:end, 1) = xi;
            M(2:stepSize:end, 2) = eta;
            M(3:stepSize:end, 3) = zeta;
            M(4:stepSize:end, 4) = xi;
            M(4:stepSize:end, 5) = eta;
            M(5:stepSize:end, 6) = eta;
            M(5:stepSize:end, 7) = zeta;
            M(6:stepSize:end, 8) = xi;
            M(6:stepSize:end, 9) = zeta;
        elseif orderAnsatzFunction == 12
            M(1:stepSize:end, 1) = xi;
            M(2:stepSize:end, 2) = eta;
            M(3:stepSize:end, 3) = zeta;
            M(4:stepSize:end, 4) = xi;
            M(4:stepSize:end, 5) = eta;
            M(5:stepSize:end, 6) = eta;
            M(5:stepSize:end, 7) = zeta;
            M(6:stepSize:end, 8) = xi;
            M(6:stepSize:end, 9) = zeta;
            M(1:stepSize:end, 10:12) = [xi .* eta; eta .* zeta; zeta .* xi].';
            M(2:stepSize:end, 10:12) = [xi .* eta; eta .* zeta; zeta .* xi].';
            M(3:stepSize:end, 10:12) = [xi .* eta; eta .* zeta; zeta .* xi].';
        else
            error('Ansatz function not implemented!');
        end
    else
        error('No ansatz functions implemented for this dimension!');
    end
end

end
