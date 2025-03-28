function [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = computeLagrangeShapeFunction(dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
%COMPUTELAGRANGESHAPEFUNCTION Computation of Lagrange shape functions.
%   This function computes Lagrange shape functions and their derivatives.

% initialize arrays
N_k_I = zeros(numberOfGausspoints, numberOfNodes);
dN_xi_k_I = zeros(dimension, numberOfGausspoints, numberOfNodes);
d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);

% allow symbolic computation of shape functions
if isa(gaussPoints, 'sym')
    N_k_I = sym(N_k_I);
    dN_xi_k_I = sym(dN_xi_k_I);
    d2N_xi_xi_k_I = sym(d2N_xi_xi_k_I);
end

% compute shape functions
if dimension == 1
    if numberOfNodes > 3
        % Automatic computation of Lagrange shape functions with more than 3 nodes
        [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = automaticComputationLagrangeShapeFunctions(dimension, numberOfNodes, gaussPoints);
    else
        % Explicit calculation of Lagrange shape functions for elements with up to 3 nodes for performance reasons
        xi = gaussPoints(1, :);
        if numberOfNodes == 1
            N_k_I(:, 1) = 1;
            dN_xi_k_I(1, :, 1) = 0;
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 2
            N_k_I(:, 1) = 1 / 2 * (1 - xi);
            N_k_I(:, 2) = 1 / 2 * (1 + xi);
            dN_xi_k_I(1, :, 1) = -1 / 2;
            dN_xi_k_I(1, :, 2) = 1 / 2;
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 3
            N_k_I(:, 1) = 1 / 2 * xi .* (xi - 1);
            N_k_I(:, 2) = 1 / 2 * xi .* (xi + 1);
            N_k_I(:, 3) = 1 - xi.^2;
            dN_xi_k_I(1, :, 1) = xi - 1 / 2;
            dN_xi_k_I(1, :, 2) = xi + 1 / 2;
            dN_xi_k_I(1, :, 3) = -2 * xi;
            d2N_xi_xi_k_I(1, 1, :, 1) = 1;
            d2N_xi_xi_k_I(1, 1, :, 2) = 1;
            d2N_xi_xi_k_I(1, 1, :, 3) = -2;
        else
            error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
        end
    end
elseif dimension == 2
    if mod(sqrt(numberOfNodes), 1) == 0 && numberOfNodes > 9
        % Automatic computation of Lagrange shape functions with more than 9 nodes
        [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = automaticComputationLagrangeShapeFunctions(dimension, numberOfNodes, gaussPoints);
    else
        % Explicit calculation of Lagrange shape functions for elements with up to 9 nodes for performance reasons
        xi = gaussPoints(1, :);
        eta = gaussPoints(2, :);
        if numberOfNodes == 1
            N_k_I(:, :) = 1;
            dN_xi_k_I(1, :, :) = 0;
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 3
            N_k_I(:, 1) = 1 - xi - eta;
            N_k_I(:, 2) = xi;
            N_k_I(:, 3) = eta;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -1;
            dN_xi_k_I(1, :, 2) = 1;
            dN_xi_k_I(1, :, 3) = 0;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = -1;
            dN_xi_k_I(2, :, 2) = 0;
            dN_xi_k_I(2, :, 3) = 1;
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 4
            N_k_I(:, 1) = (1 - xi) .* (1 - eta) / 4;
            N_k_I(:, 2) = (1 + xi) .* (1 - eta) / 4;
            N_k_I(:, 3) = (1 + xi) .* (1 + eta) / 4;
            N_k_I(:, 4) = (1 - xi) .* (1 + eta) / 4;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -(1 - eta) / 4;
            dN_xi_k_I(1, :, 2) = (1 - eta) / 4;
            dN_xi_k_I(1, :, 3) = (1 + eta) / 4;
            dN_xi_k_I(1, :, 4) = -(1 + eta) / 4;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = -(1 - xi) / 4;
            dN_xi_k_I(2, :, 2) = -(1 + xi) / 4;
            dN_xi_k_I(2, :, 3) = (1 + xi) / 4;
            dN_xi_k_I(2, :, 4) = (1 - xi) / 4;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 0;
            d2N_xi_xi_k_I(1, 1, :, 2) = 0;
            d2N_xi_xi_k_I(1, 1, :, 3) = 0;
            d2N_xi_xi_k_I(1, 1, :, 4) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = 1 / 4;
            d2N_xi_xi_k_I(1, 2, :, 2) = -1 / 4;
            d2N_xi_xi_k_I(1, 2, :, 3) = 1 / 4;
            d2N_xi_xi_k_I(1, 2, :, 4) = -1 / 4;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = 1 / 4;
            d2N_xi_xi_k_I(2, 1, :, 2) = -1 / 4;
            d2N_xi_xi_k_I(2, 1, :, 3) = 1 / 4;
            d2N_xi_xi_k_I(2, 1, :, 4) = -1 / 4;
            d2N_xi_xi_k_I(2, 2, :, 1) = 0;
            d2N_xi_xi_k_I(2, 2, :, 2) = 0;
            d2N_xi_xi_k_I(2, 2, :, 3) = 0;
            d2N_xi_xi_k_I(2, 2, :, 4) = 0;
        elseif numberOfNodes == 6
            lambda = 1 - xi - eta;
            N_k_I(:, 1) = lambda .* (2 * lambda - 1);
            N_k_I(:, 2) = xi .* (2 * xi - 1);
            N_k_I(:, 3) = eta .* (2 * eta - 1);
            N_k_I(:, 4) = 4 * xi .* lambda;
            N_k_I(:, 5) = 4 * xi .* eta;
            N_k_I(:, 6) = 4 * eta .* lambda;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = 4 * eta + 4 * xi - 3;
            dN_xi_k_I(1, :, 2) = 4 * xi - 1;
            dN_xi_k_I(1, :, 3) = 0;
            dN_xi_k_I(1, :, 4) = 4 - 8 * xi - 4 * eta;
            dN_xi_k_I(1, :, 5) = 4 * eta;
            dN_xi_k_I(1, :, 6) = -4 * eta;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = 4 * eta + 4 * xi - 3;
            dN_xi_k_I(2, :, 2) = 0;
            dN_xi_k_I(2, :, 3) = 4 * eta - 1;
            dN_xi_k_I(2, :, 4) = -4 * xi;
            dN_xi_k_I(2, :, 5) = 4 * xi;
            dN_xi_k_I(2, :, 6) = 4 - 4 * xi - 8 * eta;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 4;
            d2N_xi_xi_k_I(1, 1, :, 2) = 4;
            d2N_xi_xi_k_I(1, 1, :, 3) = 0;
            d2N_xi_xi_k_I(1, 1, :, 4) = -8;
            d2N_xi_xi_k_I(1, 1, :, 5) = 0;
            d2N_xi_xi_k_I(1, 1, :, 6) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = 4;
            d2N_xi_xi_k_I(1, 2, :, 2) = 0;
            d2N_xi_xi_k_I(1, 2, :, 3) = 0;
            d2N_xi_xi_k_I(1, 2, :, 4) = -4;
            d2N_xi_xi_k_I(1, 2, :, 5) = 4;
            d2N_xi_xi_k_I(1, 2, :, 6) = -4;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = 4;
            d2N_xi_xi_k_I(2, 1, :, 2) = 0;
            d2N_xi_xi_k_I(2, 1, :, 3) = 0;
            d2N_xi_xi_k_I(2, 1, :, 4) = -4;
            d2N_xi_xi_k_I(2, 1, :, 5) = 4;
            d2N_xi_xi_k_I(2, 1, :, 6) = -4;
            d2N_xi_xi_k_I(2, 2, :, 1) = 4;
            d2N_xi_xi_k_I(2, 2, :, 2) = 0;
            d2N_xi_xi_k_I(2, 2, :, 3) = 4;
            d2N_xi_xi_k_I(2, 2, :, 4) = 0;
            d2N_xi_xi_k_I(2, 2, :, 5) = 0;
            d2N_xi_xi_k_I(2, 2, :, 6) = -8;

        elseif numberOfNodes == 8
            N_k_I(:, 1) = (1 - xi) .* (1 - eta) .* (-xi - eta - 1) / 4;
            N_k_I(:, 2) = (1 + xi) .* (1 - eta) .* (xi - eta - 1) / 4;
            N_k_I(:, 3) = (1 + xi) .* (1 + eta) .* (xi + eta - 1) / 4;
            N_k_I(:, 4) = (1 - xi) .* (1 + eta) .* (-xi + eta - 1) / 4;
            N_k_I(:, 5) = (1 - xi.^2) .* (1 - eta) / 2;
            N_k_I(:, 6) = (1 + xi) .* (1 - eta.^2) / 2;
            N_k_I(:, 7) = (1 - xi.^2) .* (1 + eta) / 2;
            N_k_I(:, 8) = (1 - xi) .* (1 - eta.^2) / 2;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -((eta + 2 * xi) .* (eta - 1)) / 4;
            dN_xi_k_I(1, :, 2) = ((1 - eta) .* (xi - eta - 1) + (1 + xi) .* (1 - eta)) / 4;
            dN_xi_k_I(1, :, 3) = ((1 + eta) .* (xi + eta - 1) + (1 + xi) .* (1 + eta)) / 4;
            dN_xi_k_I(1, :, 4) = (-(1 + eta) .* (-xi + eta - 1) - (1 - xi) .* (1 + eta)) / 4;
            dN_xi_k_I(1, :, 5) = -xi .* (1 - eta);
            dN_xi_k_I(1, :, 6) = (1 - eta.^2) / 2;
            dN_xi_k_I(1, :, 7) = -xi .* (1 + eta);
            dN_xi_k_I(1, :, 8) = -(1 - eta.^2) / 2;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = (-(1 - xi) .* (-xi - eta - 1) - (1 - xi) .* (1 - eta)) / 4;
            dN_xi_k_I(2, :, 2) = (-(1 + xi) .* (xi - eta - 1) - (1 + xi) .* (1 - eta)) / 4;
            dN_xi_k_I(2, :, 3) = ((1 + xi) .* (xi + eta - 1) + (1 + xi) .* (1 + eta)) / 4;
            dN_xi_k_I(2, :, 4) = ((1 - xi) .* (-xi + eta - 1) + (1 - xi) .* (1 + eta)) / 4;
            dN_xi_k_I(2, :, 5) = -(1 - xi.^2) / 2;
            dN_xi_k_I(2, :, 6) = -(1 + xi) .* eta;
            dN_xi_k_I(2, :, 7) = (1 - xi.^2) / 2;
            dN_xi_k_I(2, :, 8) = -(1 - xi) .* eta;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 0.5 .* (1 - eta);
            d2N_xi_xi_k_I(1, 1, :, 2) = 0.5 .* (1 - eta);
            d2N_xi_xi_k_I(1, 1, :, 3) = 0.5 .* (1 + eta);
            d2N_xi_xi_k_I(1, 1, :, 4) = 0.5 .* (1 + eta);
            d2N_xi_xi_k_I(1, 1, :, 5) = eta - 1;
            d2N_xi_xi_k_I(1, 1, :, 6) = 0;
            d2N_xi_xi_k_I(1, 1, :, 7) = -eta - 1;
            d2N_xi_xi_k_I(1, 1, :, 8) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = 1 / 4 - xi / 2 - eta / 2;
            d2N_xi_xi_k_I(1, 2, :, 2) = -0.5 .* xi + 0.5 .* eta - 0.25;
            d2N_xi_xi_k_I(1, 2, :, 3) = 0.5 .* xi + 0.5 .* eta + 0.25;
            d2N_xi_xi_k_I(1, 2, :, 4) = 0.5 .* xi - 0.5 .* eta - 0.25;
            d2N_xi_xi_k_I(1, 2, :, 5) = xi;
            d2N_xi_xi_k_I(1, 2, :, 6) = -eta;
            d2N_xi_xi_k_I(1, 2, :, 7) = -xi;
            d2N_xi_xi_k_I(1, 2, :, 8) = eta;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = -0.5 .* xi - 0.5 .* eta + 0.25;
            d2N_xi_xi_k_I(2, 1, :, 2) = -0.5 .* xi + 0.5 .* eta - 0.25;
            d2N_xi_xi_k_I(2, 1, :, 3) = 0.5 .* xi + 0.5 .* eta + 0.25;
            d2N_xi_xi_k_I(2, 1, :, 4) = 0.5 .* xi - 0.5 .* eta - 0.25;
            d2N_xi_xi_k_I(2, 1, :, 5) = xi;
            d2N_xi_xi_k_I(2, 1, :, 6) = -eta;
            d2N_xi_xi_k_I(2, 1, :, 7) = -xi;
            d2N_xi_xi_k_I(2, 1, :, 8) = eta;
            d2N_xi_xi_k_I(2, 2, :, 1) = 0.5 .* (1 - xi);
            d2N_xi_xi_k_I(2, 2, :, 2) = 0.5 .* (1 + xi);
            d2N_xi_xi_k_I(2, 2, :, 3) = 0.5 .* (1 + xi);
            d2N_xi_xi_k_I(2, 2, :, 4) = 0.5 .* (1 - xi);
            d2N_xi_xi_k_I(2, 2, :, 5) = 0;
            d2N_xi_xi_k_I(2, 2, :, 6) = -xi - 1;
            d2N_xi_xi_k_I(2, 2, :, 7) = 0;
            d2N_xi_xi_k_I(2, 2, :, 8) = xi - 1;

        elseif numberOfNodes == 9
            N_k_I(:, 1) = 1 / 4 * (xi.^2 - xi) .* (eta.^2 - eta);
            N_k_I(:, 2) = 1 / 4 * (xi.^2 + xi) .* (eta.^2 - eta);
            N_k_I(:, 3) = 1 / 4 * (xi.^2 + xi) .* (eta.^2 + eta);
            N_k_I(:, 4) = 1 / 4 * (xi.^2 - xi) .* (eta.^2 + eta);
            N_k_I(:, 5) = 1 / 2 * (1 - xi.^2) .* eta .* (eta - 1);
            N_k_I(:, 6) = 1 / 2 * xi .* (xi + 1) .* (1 - eta.^2);
            N_k_I(:, 7) = 1 / 2 * (1 - xi.^2) .* eta .* (eta + 1);
            N_k_I(:, 8) = 1 / 2 * xi .* (xi - 1) .* (1 - eta.^2);
            N_k_I(:, 9) = (1 - xi.^2) .* (1 - eta.^2);
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = 1 / 4 * (2 * xi - 1) .* (eta.^2 - eta);
            dN_xi_k_I(1, :, 2) = 1 / 4 * (2 * xi + 1) .* (eta.^2 - eta);
            dN_xi_k_I(1, :, 3) = 1 / 4 * (2 * xi + 1) .* (eta.^2 + eta);
            dN_xi_k_I(1, :, 4) = 1 / 4 * (2 * xi - 1) .* (eta.^2 + eta);
            dN_xi_k_I(1, :, 5) = -xi .* eta .* (eta - 1);
            dN_xi_k_I(1, :, 6) = 1 / 2 * (2 * xi + 1) .* (1 - eta.^2);
            dN_xi_k_I(1, :, 7) = -xi .* eta .* (eta + 1);
            dN_xi_k_I(1, :, 8) = 1 / 2 * (2 * xi - 1) .* (1 - eta.^2);
            dN_xi_k_I(1, :, 9) = -2 * xi .* (1 - eta.^2);
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = 1 / 4 * (xi.^2 - xi) .* (2 * eta - 1);
            dN_xi_k_I(2, :, 2) = 1 / 4 * (xi.^2 + xi) .* (2 * eta - 1);
            dN_xi_k_I(2, :, 3) = 1 / 4 * (xi.^2 + xi) .* (2 * eta + 1);
            dN_xi_k_I(2, :, 4) = 1 / 4 * (xi.^2 - xi) .* (2 * eta + 1);
            dN_xi_k_I(2, :, 5) = 1 / 2 * (1 - xi.^2) .* (2 * eta - 1);
            dN_xi_k_I(2, :, 6) = -xi .* (xi + 1) .* eta;
            dN_xi_k_I(2, :, 7) = 1 / 2 * (1 - xi.^2) .* (2 * eta + 1);
            dN_xi_k_I(2, :, 8) = -xi .* (xi - 1) .* eta;
            dN_xi_k_I(2, :, 9) = -2 * (1 - xi.^2) .* eta;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = -2 * eta .* (1 - eta) / 4;
            d2N_xi_xi_k_I(1, 1, :, 2) = -2 * eta .* (1 - eta) / 4;
            d2N_xi_xi_k_I(1, 1, :, 3) = 2 * eta .* (1 + eta) / 4;
            d2N_xi_xi_k_I(1, 1, :, 4) = 2 * eta .* (1 + eta) / 4;
            d2N_xi_xi_k_I(1, 1, :, 5) = 2 * eta .* (1 - eta) / 2;
            d2N_xi_xi_k_I(1, 1, :, 6) = 2 * (1 - eta.^2) / 2;
            d2N_xi_xi_k_I(1, 1, :, 7) = -2 * eta .* (1 + eta) / 2;
            d2N_xi_xi_k_I(1, 1, :, 8) = 2 * (1 - eta.^2) / 2;
            d2N_xi_xi_k_I(1, 1, :, 9) = -2 * (1 - eta.^2);
            d2N_xi_xi_k_I(1, 2, :, 1) = (1 - 2 * xi) .* (1 - 2 * eta) / 4;
            d2N_xi_xi_k_I(1, 2, :, 2) = -(1 + 2 * xi) .* (1 - 2 * eta) / 4;
            d2N_xi_xi_k_I(1, 2, :, 3) = (1 + 2 * xi) .* (1 + 2 * eta) / 4;
            d2N_xi_xi_k_I(1, 2, :, 4) = -(1 - 2 * xi) .* (1 + 2 * eta) / 4;
            d2N_xi_xi_k_I(1, 2, :, 5) = 2 * xi .* (1 - 2 * eta) / 2;
            d2N_xi_xi_k_I(1, 2, :, 6) = (1 + 2 * xi) .* (-2 * eta) / 2;
            d2N_xi_xi_k_I(1, 2, :, 7) = -2 * xi .* (1 + 2 * eta) / 2;
            d2N_xi_xi_k_I(1, 2, :, 8) = -(1 - 2 * xi) .* (-2 * eta) / 2;
            d2N_xi_xi_k_I(1, 2, :, 9) = -2 * xi .* (-2 * eta);
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = (1 - 2 * xi) .* (1 - 2 * eta) / 4;
            d2N_xi_xi_k_I(2, 1, :, 2) = -(1 + 2 * xi) .* (1 - 2 * eta) / 4;
            d2N_xi_xi_k_I(2, 1, :, 3) = (1 + 2 * xi) .* (1 + 2 * eta) / 4;
            d2N_xi_xi_k_I(2, 1, :, 4) = -(1 - 2 * xi) .* (1 + 2 * eta) / 4;
            d2N_xi_xi_k_I(2, 1, :, 5) = 2 * xi .* (1 - 2 * eta) / 2;
            d2N_xi_xi_k_I(2, 1, :, 6) = (1 + 2 * xi) .* (-2 * eta) / 2;
            d2N_xi_xi_k_I(2, 1, :, 7) = -2 * xi .* (1 + 2 * eta) / 2;
            d2N_xi_xi_k_I(2, 1, :, 8) = -(1 - 2 * xi) .* (-2 * eta) / 2;
            d2N_xi_xi_k_I(2, 1, :, 9) = -2 * xi .* (-2 * eta);
            d2N_xi_xi_k_I(2, 2, :, 1) = -2 * xi .* (1 - xi) / 4;
            d2N_xi_xi_k_I(2, 2, :, 2) = 2 * xi .* (1 + xi) / 4;
            d2N_xi_xi_k_I(2, 2, :, 3) = 2 * xi .* (1 + xi) / 4;
            d2N_xi_xi_k_I(2, 2, :, 4) = -2 * xi .* (1 - xi) / 4;
            d2N_xi_xi_k_I(2, 2, :, 5) = 2 * (1 - xi.^2) / 2;
            d2N_xi_xi_k_I(2, 2, :, 6) = -2 * xi .* (1 + xi) / 2;
            d2N_xi_xi_k_I(2, 2, :, 7) = 2 * (1 - xi.^2) / 2;
            d2N_xi_xi_k_I(2, 2, :, 8) = 2 * xi .* (1 - xi) / 2;
            d2N_xi_xi_k_I(2, 2, :, 9) = -2 * (1 - xi.^2);
        else
            error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
        end
    end
elseif dimension == 3
    if mod(nthroot(numberOfNodes, dimension), 1) == 0 && numberOfNodes > 27
        % Automatic computation of Lagrange shape functions with more than 27 nodes
        [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = automaticComputationLagrangeShapeFunctions(dimension, numberOfNodes, gaussPoints);
    else
        % Explicit calculation of Lagrange shape functions for elements with up to 27 nodes for performance reasons
        xi = gaussPoints(1, :);
        eta = gaussPoints(2, :);
        zeta = gaussPoints(3, :);
        if numberOfNodes == 1
            N_k_I(1:numberOfGausspoints, 1) = 1;
            dN_xi_k_I(1, 1:numberOfGausspoints, 1) = 0;
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 4
            N_k_I(:, 1) = 1 - xi - eta - zeta;
            N_k_I(:, 2) = xi;
            N_k_I(:, 3) = eta;
            N_k_I(:, 4) = zeta;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -1;
            dN_xi_k_I(1, :, 2) = 1;
            dN_xi_k_I(1, :, 3) = 0;
            dN_xi_k_I(1, :, 4) = 0;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = -1;
            dN_xi_k_I(2, :, 2) = 0;
            dN_xi_k_I(2, :, 3) = 1;
            dN_xi_k_I(2, :, 4) = 0;
            % derivative wrt zeta
            dN_xi_k_I(3, :, 1) = -1;
            dN_xi_k_I(3, :, 2) = 0;
            dN_xi_k_I(3, :, 3) = 0;
            dN_xi_k_I(3, :, 4) = 1;
            % second derivative
            d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);
        elseif numberOfNodes == 8
            N_k_I(:, 1) = (1 - xi) .* (1 - eta) .* (1 - zeta) / 8;
            N_k_I(:, 2) = (1 + xi) .* (1 - eta) .* (1 - zeta) / 8;
            N_k_I(:, 3) = (1 + xi) .* (1 + eta) .* (1 - zeta) / 8;
            N_k_I(:, 4) = (1 - xi) .* (1 + eta) .* (1 - zeta) / 8;
            N_k_I(:, 5) = (1 - xi) .* (1 - eta) .* (1 + zeta) / 8;
            N_k_I(:, 6) = (1 + xi) .* (1 - eta) .* (1 + zeta) / 8;
            N_k_I(:, 7) = (1 + xi) .* (1 + eta) .* (1 + zeta) / 8;
            N_k_I(:, 8) = (1 - xi) .* (1 + eta) .* (1 + zeta) / 8;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -(1 - eta) .* (1 - zeta) / 8;
            dN_xi_k_I(1, :, 2) = (1 - eta) .* (1 - zeta) / 8;
            dN_xi_k_I(1, :, 3) = (1 + eta) .* (1 - zeta) / 8;
            dN_xi_k_I(1, :, 4) = -(1 + eta) .* (1 - zeta) / 8;
            dN_xi_k_I(1, :, 5) = -(1 - eta) .* (1 + zeta) / 8;
            dN_xi_k_I(1, :, 6) = (1 - eta) .* (1 + zeta) / 8;
            dN_xi_k_I(1, :, 7) = (1 + eta) .* (1 + zeta) / 8;
            dN_xi_k_I(1, :, 8) = -(1 + eta) .* (1 + zeta) / 8;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = -(1 - xi) .* (1 - zeta) / 8;
            dN_xi_k_I(2, :, 2) = -(1 + xi) .* (1 - zeta) / 8;
            dN_xi_k_I(2, :, 3) = (1 + xi) .* (1 - zeta) / 8;
            dN_xi_k_I(2, :, 4) = (1 - xi) .* (1 - zeta) / 8;
            dN_xi_k_I(2, :, 5) = -(1 - xi) .* (1 + zeta) / 8;
            dN_xi_k_I(2, :, 6) = -(1 + xi) .* (1 + zeta) / 8;
            dN_xi_k_I(2, :, 7) = (1 + xi) .* (1 + zeta) / 8;
            dN_xi_k_I(2, :, 8) = (1 - xi) .* (1 + zeta) / 8;
            % derivative wrt zeta
            dN_xi_k_I(3, :, 1) = -(1 - xi) .* (1 - eta) / 8;
            dN_xi_k_I(3, :, 2) = -(1 + xi) .* (1 - eta) / 8;
            dN_xi_k_I(3, :, 3) = -(1 + xi) .* (1 + eta) / 8;
            dN_xi_k_I(3, :, 4) = -(1 - xi) .* (1 + eta) / 8;
            dN_xi_k_I(3, :, 5) = (1 - xi) .* (1 - eta) / 8;
            dN_xi_k_I(3, :, 6) = (1 + xi) .* (1 - eta) / 8;
            dN_xi_k_I(3, :, 7) = (1 + xi) .* (1 + eta) / 8;
            dN_xi_k_I(3, :, 8) = (1 - xi) .* (1 + eta) / 8;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 0;
            d2N_xi_xi_k_I(1, 1, :, 2) = 0;
            d2N_xi_xi_k_I(1, 1, :, 3) = 0;
            d2N_xi_xi_k_I(1, 1, :, 4) = 0;
            d2N_xi_xi_k_I(1, 1, :, 5) = 0;
            d2N_xi_xi_k_I(1, 1, :, 6) = 0;
            d2N_xi_xi_k_I(1, 1, :, 7) = 0;
            d2N_xi_xi_k_I(1, 1, :, 8) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = (1 - zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 2) = -(1 - zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 3) = (1 - zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 4) = -(1 - zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 5) = (1 + zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 6) = -(1 + zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 7) = (1 + zeta) / 8;
            d2N_xi_xi_k_I(1, 2, :, 8) = -(1 + zeta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 1) = (1 - eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 2) = -(1 - eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 3) = -(1 + eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 4) = (1 + eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 5) = -(1 - eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 6) = (1 - eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 7) = (1 + eta) / 8;
            d2N_xi_xi_k_I(1, 3, :, 8) = -(1 + eta) / 8;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = (1 - zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 2) = -(1 - zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 3) = (1 - zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 4) = -(1 - zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 5) = (1 + zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 6) = -(1 + zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 7) = (1 + zeta) / 8;
            d2N_xi_xi_k_I(2, 1, :, 8) = -(1 + zeta) / 8;
            d2N_xi_xi_k_I(2, 2, :, 1) = 0;
            d2N_xi_xi_k_I(2, 2, :, 2) = 0;
            d2N_xi_xi_k_I(2, 2, :, 3) = 0;
            d2N_xi_xi_k_I(2, 2, :, 4) = 0;
            d2N_xi_xi_k_I(2, 2, :, 5) = 0;
            d2N_xi_xi_k_I(2, 2, :, 6) = 0;
            d2N_xi_xi_k_I(2, 2, :, 7) = 0;
            d2N_xi_xi_k_I(2, 2, :, 8) = 0;
            d2N_xi_xi_k_I(2, 3, :, 1) = (1 - xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 2) = (1 + xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 3) = -(1 + xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 4) = -(1 - xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 5) = -(1 - xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 6) = -(1 + xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 7) = (1 + xi) / 8;
            d2N_xi_xi_k_I(2, 3, :, 8) = (1 - xi) / 8;
            % second derivative wrt zeta
            d2N_xi_xi_k_I(3, 1, :, 1) = (1 - eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 2) = -(1 - eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 3) = -(1 + eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 4) = (1 + eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 5) = -(1 - eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 6) = (1 - eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 7) = (1 + eta) / 8;
            d2N_xi_xi_k_I(3, 1, :, 8) = -(1 + eta) / 8;
            d2N_xi_xi_k_I(3, 2, :, 1) = (1 - xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 2) = (1 + xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 3) = -(1 + xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 4) = -(1 - xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 5) = -(1 - xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 6) = -(1 + xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 7) = (1 + xi) / 8;
            d2N_xi_xi_k_I(3, 2, :, 8) = (1 - xi) / 8;
            d2N_xi_xi_k_I(3, 3, :, 1) = 0;
            d2N_xi_xi_k_I(3, 3, :, 2) = 0;
            d2N_xi_xi_k_I(3, 3, :, 3) = 0;
            d2N_xi_xi_k_I(3, 3, :, 4) = 0;
            d2N_xi_xi_k_I(3, 3, :, 5) = 0;
            d2N_xi_xi_k_I(3, 3, :, 6) = 0;
            d2N_xi_xi_k_I(3, 3, :, 7) = 0;
            d2N_xi_xi_k_I(3, 3, :, 8) = 0;
        elseif numberOfNodes == 10
            lambda = 1 - xi - eta - zeta;
            N_k_I(:, 1) = lambda .* (2 * lambda - 1);
            N_k_I(:, 2) = xi .* (2 * xi - 1);
            N_k_I(:, 3) = eta .* (2 * eta - 1);
            N_k_I(:, 4) = zeta .* (2 * zeta - 1);
            N_k_I(:, 5) = 4 * xi .* lambda;
            N_k_I(:, 6) = 4 * xi .* eta;
            N_k_I(:, 7) = 4 * eta .* lambda;
            N_k_I(:, 8) = 4 * zeta .* lambda;
            N_k_I(:, 9) = 4 * xi .* zeta;
            N_k_I(:, 10) = 4 * eta .* zeta;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = -3 + 4 * xi + 4 * eta + 4 * zeta;
            dN_xi_k_I(1, :, 2) = 4 * xi - 1;
            dN_xi_k_I(1, :, 3) = 0;
            dN_xi_k_I(1, :, 4) = 0;
            dN_xi_k_I(1, :, 5) = 4 - 8 * xi - 4 * eta - 4 * zeta;
            dN_xi_k_I(1, :, 6) = 4 * eta;
            dN_xi_k_I(1, :, 7) = -4 * eta;
            dN_xi_k_I(1, :, 8) = -4 * zeta;
            dN_xi_k_I(1, :, 9) = 4 * zeta;
            dN_xi_k_I(1, :, 10) = 0;
            % Derivative wrt eta
            dN_xi_k_I(2, :, 1) = -3 + 4 * xi + 4 * eta + 4 * zeta;
            dN_xi_k_I(2, :, 2) = 0;
            dN_xi_k_I(2, :, 3) = 4 * eta - 1;
            dN_xi_k_I(2, :, 4) = 0;
            dN_xi_k_I(2, :, 5) = -4 * xi;
            dN_xi_k_I(2, :, 6) = 4 * xi;
            dN_xi_k_I(2, :, 7) = 4 - 4 * xi - 8 * eta - 4 * zeta;
            dN_xi_k_I(2, :, 8) = -4 * zeta;
            dN_xi_k_I(2, :, 9) = 0;
            dN_xi_k_I(2, :, 10) = 4 * zeta;
            % derivative wrt zeta
            dN_xi_k_I(3, :, 1) = -3 + 4 * xi + 4 * eta + 4 * zeta;
            dN_xi_k_I(3, :, 2) = 0;
            dN_xi_k_I(3, :, 3) = 0;
            dN_xi_k_I(3, :, 4) = 4 * zeta - 1;
            dN_xi_k_I(3, :, 5) = -4 * xi;
            dN_xi_k_I(3, :, 6) = 0;
            dN_xi_k_I(3, :, 7) = -4 * eta;
            dN_xi_k_I(3, :, 8) = 4 - 4 * xi - 4 * eta - 8 * zeta;
            dN_xi_k_I(3, :, 9) = 4 * xi;
            dN_xi_k_I(3, :, 10) = 4 * eta;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 4;
            d2N_xi_xi_k_I(1, 1, :, 2) = 4;
            d2N_xi_xi_k_I(1, 1, :, 3) = 0;
            d2N_xi_xi_k_I(1, 1, :, 4) = 0;
            d2N_xi_xi_k_I(1, 1, :, 5) = -8;
            d2N_xi_xi_k_I(1, 1, :, 6) = 0;
            d2N_xi_xi_k_I(1, 1, :, 7) = 0;
            d2N_xi_xi_k_I(1, 1, :, 8) = 0;
            d2N_xi_xi_k_I(1, 1, :, 9) = 0;
            d2N_xi_xi_k_I(1, 1, :, 10) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = 4;
            d2N_xi_xi_k_I(1, 2, :, 2) = 0;
            d2N_xi_xi_k_I(1, 2, :, 3) = 0;
            d2N_xi_xi_k_I(1, 2, :, 4) = 0;
            d2N_xi_xi_k_I(1, 2, :, 5) = -4;
            d2N_xi_xi_k_I(1, 2, :, 6) = 4;
            d2N_xi_xi_k_I(1, 2, :, 7) = -4;
            d2N_xi_xi_k_I(1, 2, :, 8) = 0;
            d2N_xi_xi_k_I(1, 2, :, 9) = 0;
            d2N_xi_xi_k_I(1, 2, :, 10) = 0;
            d2N_xi_xi_k_I(1, 3, :, 1) = 4;
            d2N_xi_xi_k_I(1, 3, :, 2) = 0;
            d2N_xi_xi_k_I(1, 3, :, 3) = 0;
            d2N_xi_xi_k_I(1, 3, :, 4) = 0;
            d2N_xi_xi_k_I(1, 3, :, 5) = -4;
            d2N_xi_xi_k_I(1, 3, :, 6) = 0;
            d2N_xi_xi_k_I(1, 3, :, 7) = 0;
            d2N_xi_xi_k_I(1, 3, :, 8) = -4;
            d2N_xi_xi_k_I(1, 3, :, 9) = 4;
            d2N_xi_xi_k_I(1, 3, :, 10) = 0;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = 4;
            d2N_xi_xi_k_I(2, 1, :, 2) = 0;
            d2N_xi_xi_k_I(2, 1, :, 3) = 0;
            d2N_xi_xi_k_I(2, 1, :, 4) = 0;
            d2N_xi_xi_k_I(2, 1, :, 5) = -4;
            d2N_xi_xi_k_I(2, 1, :, 6) = 4;
            d2N_xi_xi_k_I(2, 1, :, 7) = -4;
            d2N_xi_xi_k_I(2, 1, :, 8) = 0;
            d2N_xi_xi_k_I(2, 1, :, 9) = 0;
            d2N_xi_xi_k_I(2, 1, :, 10) = 0;
            d2N_xi_xi_k_I(2, 2, :, 1) = 4;
            d2N_xi_xi_k_I(2, 2, :, 2) = 0;
            d2N_xi_xi_k_I(2, 2, :, 3) = 4;
            d2N_xi_xi_k_I(2, 2, :, 4) = 0;
            d2N_xi_xi_k_I(2, 2, :, 5) = 0;
            d2N_xi_xi_k_I(2, 2, :, 6) = 0;
            d2N_xi_xi_k_I(2, 2, :, 7) = -8;
            d2N_xi_xi_k_I(2, 2, :, 8) = 0;
            d2N_xi_xi_k_I(2, 2, :, 9) = 0;
            d2N_xi_xi_k_I(2, 2, :, 10) = 0;
            d2N_xi_xi_k_I(2, 3, :, 1) = 4;
            d2N_xi_xi_k_I(2, 3, :, 2) = 0;
            d2N_xi_xi_k_I(2, 3, :, 3) = 0;
            d2N_xi_xi_k_I(2, 3, :, 4) = 0;
            d2N_xi_xi_k_I(2, 3, :, 5) = 0;
            d2N_xi_xi_k_I(2, 3, :, 6) = 0;
            d2N_xi_xi_k_I(2, 3, :, 7) = -4;
            d2N_xi_xi_k_I(2, 3, :, 8) = -4;
            d2N_xi_xi_k_I(2, 3, :, 9) = 0;
            d2N_xi_xi_k_I(2, 3, :, 10) = 4;
            % second derivative wrt zeta
            d2N_xi_xi_k_I(3, 1, :, 1) = 4;
            d2N_xi_xi_k_I(3, 1, :, 2) = 0;
            d2N_xi_xi_k_I(3, 1, :, 3) = 0;
            d2N_xi_xi_k_I(3, 1, :, 4) = 0;
            d2N_xi_xi_k_I(3, 1, :, 5) = -4;
            d2N_xi_xi_k_I(3, 1, :, 6) = 0;
            d2N_xi_xi_k_I(3, 1, :, 7) = 0;
            d2N_xi_xi_k_I(3, 1, :, 8) = -4;
            d2N_xi_xi_k_I(3, 1, :, 9) = 4;
            d2N_xi_xi_k_I(3, 1, :, 10) = 0;
            d2N_xi_xi_k_I(3, 2, :, 1) = 4;
            d2N_xi_xi_k_I(3, 2, :, 2) = 0;
            d2N_xi_xi_k_I(3, 2, :, 3) = 0;
            d2N_xi_xi_k_I(3, 2, :, 4) = 0;
            d2N_xi_xi_k_I(3, 2, :, 5) = 0;
            d2N_xi_xi_k_I(3, 2, :, 6) = 0;
            d2N_xi_xi_k_I(3, 2, :, 7) = -4;
            d2N_xi_xi_k_I(3, 2, :, 8) = -4;
            d2N_xi_xi_k_I(3, 2, :, 9) = 0;
            d2N_xi_xi_k_I(3, 2, :, 10) = 4;
            d2N_xi_xi_k_I(3, 3, :, 1) = 4;
            d2N_xi_xi_k_I(3, 3, :, 2) = 0;
            d2N_xi_xi_k_I(3, 3, :, 3) = 0;
            d2N_xi_xi_k_I(3, 3, :, 4) = 4;
            d2N_xi_xi_k_I(3, 3, :, 5) = 0;
            d2N_xi_xi_k_I(3, 3, :, 6) = 0;
            d2N_xi_xi_k_I(3, 3, :, 7) = 0;
            d2N_xi_xi_k_I(3, 3, :, 8) = -8;
            d2N_xi_xi_k_I(3, 3, :, 9) = 0;
            d2N_xi_xi_k_I(3, 3, :, 10) = 0;
        elseif numberOfNodes == 20
            N_k_I(:, 1) = (1 - xi) .* (1 - eta) .* (1 - zeta) .* (-xi - eta - zeta - 2) / 8;
            N_k_I(:, 2) = (1 + xi) .* (1 - eta) .* (1 - zeta) .* (xi - eta - zeta - 2) / 8;
            N_k_I(:, 3) = (1 + xi) .* (1 + eta) .* (1 - zeta) .* (xi + eta - zeta - 2) / 8;
            N_k_I(:, 4) = (1 - xi) .* (1 + eta) .* (1 - zeta) .* (-xi + eta - zeta - 2) / 8;
            N_k_I(:, 5) = (1 - xi) .* (1 - eta) .* (1 + zeta) .* (-xi - eta + zeta - 2) / 8;
            N_k_I(:, 6) = (1 + xi) .* (1 - eta) .* (1 + zeta) .* (xi - eta + zeta - 2) / 8;
            N_k_I(:, 7) = (1 + xi) .* (1 + eta) .* (1 + zeta) .* (xi + eta + zeta - 2) / 8;
            N_k_I(:, 8) = (1 - xi) .* (1 + eta) .* (1 + zeta) .* (-xi + eta + zeta - 2) / 8;
            N_k_I(:, 9) = (1 - xi.^2) .* (1 - eta) .* (1 - zeta) / 4;
            N_k_I(:, 10) = (1 + xi) .* (1 - eta.^2) .* (1 - zeta) / 4;
            N_k_I(:, 11) = (1 - xi.^2) .* (1 + eta) .* (1 - zeta) / 4;
            N_k_I(:, 12) = (1 - xi) .* (1 - eta.^2) .* (1 - zeta) / 4;
            N_k_I(:, 13) = (1 - xi.^2) .* (1 - eta) .* (1 + zeta) / 4;
            N_k_I(:, 14) = (1 + xi) .* (1 - eta.^2) .* (1 + zeta) / 4;
            N_k_I(:, 15) = (1 - xi.^2) .* (1 + eta) .* (1 + zeta) / 4;
            N_k_I(:, 16) = (1 - xi) .* (1 - eta.^2) .* (1 + zeta) / 4;
            N_k_I(:, 17) = (1 - xi) .* (1 - eta) .* (1 - zeta.^2) / 4;
            N_k_I(:, 18) = (1 + xi) .* (1 - eta) .* (1 - zeta.^2) / 4;
            N_k_I(:, 19) = (1 + xi) .* (1 + eta) .* (1 - zeta.^2) / 4;
            N_k_I(:, 20) = (1 - xi) .* (1 + eta) .* (1 - zeta.^2) / 4;
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = (-(1 - eta) .* (1 - zeta) .* (-xi - eta - zeta - 2) - (1 - xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(1, :, 2) = ((1 - eta) .* (1 - zeta) .* (xi - eta - zeta - 2) + (1 + xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(1, :, 3) = ((1 + eta) .* (1 - zeta) .* (xi + eta - zeta - 2) + (1 + xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(1, :, 4) = (-(1 + eta) .* (1 - zeta) .* (-xi + eta - zeta - 2) - (1 - xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(1, :, 5) = (-(1 - eta) .* (1 + zeta) .* (-xi - eta + zeta - 2) - (1 - xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(1, :, 6) = ((1 - eta) .* (1 + zeta) .* (xi - eta + zeta - 2) + (1 + xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(1, :, 7) = ((1 + eta) .* (1 + zeta) .* (xi + eta + zeta - 2) + (1 + xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(1, :, 8) = (-(1 + eta) .* (1 + zeta) .* (-xi + eta + zeta - 2) - (1 - xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(1, :, 9) = -xi .* (1 - eta) .* (1 - zeta) / 2;
            dN_xi_k_I(1, :, 10) = (1 - eta.^2) .* (1 - zeta) / 4;
            dN_xi_k_I(1, :, 11) = -xi .* (1 + eta) .* (1 - zeta) / 2;
            dN_xi_k_I(1, :, 12) = -(1 - eta.^2) .* (1 - zeta) / 4;
            dN_xi_k_I(1, :, 13) = -xi .* (1 - eta) .* (1 + zeta) / 2;
            dN_xi_k_I(1, :, 14) = (1 - eta.^2) .* (1 + zeta) / 4;
            dN_xi_k_I(1, :, 15) = -xi .* (1 + eta) .* (1 + zeta) / 2;
            dN_xi_k_I(1, :, 16) = -(1 - eta.^2) .* (1 + zeta) / 4;
            dN_xi_k_I(1, :, 17) = -(1 - eta) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(1, :, 18) = (1 - eta) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(1, :, 19) = (1 + eta) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(1, :, 20) = -(1 + eta) .* (1 - zeta.^2) / 4;
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = (-(1 - xi) .* (1 - zeta) .* (-xi - eta - zeta - 2) - (1 - xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(2, :, 2) = (-(1 + xi) .* (1 - zeta) .* (xi - eta - zeta - 2) - (1 + xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(2, :, 3) = ((1 + xi) .* (1 - zeta) .* (xi + eta - zeta - 2) + (1 + xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(2, :, 4) = ((1 - xi) .* (1 - zeta) .* (-xi + eta - zeta - 2) + (1 - xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(2, :, 5) = (-(1 - xi) .* (1 + zeta) .* (-xi - eta + zeta - 2) - (1 - xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(2, :, 6) = (-(1 + xi) .* (1 + zeta) .* (xi - eta + zeta - 2) - (1 + xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(2, :, 7) = ((1 + xi) .* (1 + zeta) .* (xi + eta + zeta - 2) + (1 + xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(2, :, 8) = ((1 - xi) .* (1 + zeta) .* (-xi + eta + zeta - 2) + (1 - xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(2, :, 9) = -(1 - xi.^2) .* (1 - zeta) / 4;
            dN_xi_k_I(2, :, 10) = -(1 + xi) .* eta .* (1 - zeta) / 2;
            dN_xi_k_I(2, :, 11) = (1 - xi.^2) .* (1 - zeta) / 4;
            dN_xi_k_I(2, :, 12) = -(1 - xi) .* eta .* (1 - zeta) / 2;
            dN_xi_k_I(2, :, 13) = -(1 - xi.^2) .* (1 + zeta) / 4;
            dN_xi_k_I(2, :, 14) = -(1 + xi) .* eta .* (1 + zeta) / 2;
            dN_xi_k_I(2, :, 15) = (1 - xi.^2) .* (1 + zeta) / 4;
            dN_xi_k_I(2, :, 16) = -(1 - xi) .* eta .* (1 + zeta) / 2;
            dN_xi_k_I(2, :, 17) = -(1 - xi) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(2, :, 18) = -(1 + xi) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(2, :, 19) = (1 + xi) .* (1 - zeta.^2) / 4;
            dN_xi_k_I(2, :, 20) = (1 - xi) .* (1 - zeta.^2) / 4;
            % derivative wrt zeta
            dN_xi_k_I(3, :, 1) = (-(1 - xi) .* (1 - eta) .* (-xi - eta - zeta - 2) - (1 - xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(3, :, 2) = (-(1 + xi) .* (1 - eta) .* (xi - eta - zeta - 2) - (1 + xi) .* (1 - eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(3, :, 3) = (-(1 + xi) .* (1 + eta) .* (xi + eta - zeta - 2) - (1 + xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(3, :, 4) = (-(1 - xi) .* (1 + eta) .* (-xi + eta - zeta - 2) - (1 - xi) .* (1 + eta) .* (1 - zeta)) / 8;
            dN_xi_k_I(3, :, 5) = ((1 - xi) .* (1 - eta) .* (-xi - eta + zeta - 2) + (1 - xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(3, :, 6) = ((1 + xi) .* (1 - eta) .* (xi - eta + zeta - 2) + (1 + xi) .* (1 - eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(3, :, 7) = ((1 + xi) .* (1 + eta) .* (xi + eta + zeta - 2) + (1 + xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(3, :, 8) = ((1 - xi) .* (1 + eta) .* (-xi + eta + zeta - 2) + (1 - xi) .* (1 + eta) .* (1 + zeta)) / 8;
            dN_xi_k_I(3, :, 9) = -(1 - xi.^2) .* (1 - eta) / 4;
            dN_xi_k_I(3, :, 10) = -(1 + xi) .* (1 - eta.^2) / 4;
            dN_xi_k_I(3, :, 11) = -(1 - xi.^2) .* (1 + eta) / 4;
            dN_xi_k_I(3, :, 12) = -(1 - xi) .* (1 - eta.^2) / 4;
            dN_xi_k_I(3, :, 13) = (1 - xi.^2) .* (1 - eta) / 4;
            dN_xi_k_I(3, :, 14) = (1 + xi) .* (1 - eta.^2) / 4;
            dN_xi_k_I(3, :, 15) = (1 - xi.^2) .* (1 + eta) / 4;
            dN_xi_k_I(3, :, 16) = (1 - xi) .* (1 - eta.^2) / 4;
            dN_xi_k_I(3, :, 17) = -(1 - xi) .* (1 - eta) .* zeta / 2;
            dN_xi_k_I(3, :, 18) = -(1 + xi) .* (1 - eta) .* zeta / 2;
            dN_xi_k_I(3, :, 19) = -(1 + xi) .* (1 + eta) .* zeta / 2;
            dN_xi_k_I(3, :, 20) = -(1 - xi) .* (1 + eta) .* zeta / 2;
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = 0.25 * (eta .* zeta - eta - zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 2) = ((eta - 1) .* (zeta - 1)) / 4;
            d2N_xi_xi_k_I(1, 1, :, 3) = 0.25 * (-eta .* zeta + eta - zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 4) = -((eta + 1) .* (zeta - 1)) / 4;
            d2N_xi_xi_k_I(1, 1, :, 5) = 0.25 * (-eta .* zeta - eta + zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 6) = 0.25 * (-eta .* zeta - eta + zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 7) = 0.25 * (eta .* zeta + eta + zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 8) = 0.25 * (eta .* zeta + eta + zeta + 1);
            d2N_xi_xi_k_I(1, 1, :, 9) = 0.5 * (-eta .* zeta + eta + zeta - 1);
            d2N_xi_xi_k_I(1, 1, :, 10) = 0;
            d2N_xi_xi_k_I(1, 1, :, 11) = 0.5 * (eta .* zeta - eta + zeta - 1);
            d2N_xi_xi_k_I(1, 1, :, 12) = 0;
            d2N_xi_xi_k_I(1, 1, :, 13) = ((eta - 1) .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(1, 1, :, 14) = 0;
            d2N_xi_xi_k_I(1, 1, :, 15) = 0.5 * (-eta .* zeta - eta - zeta - 1);
            d2N_xi_xi_k_I(1, 1, :, 16) = 0;
            d2N_xi_xi_k_I(1, 1, :, 17) = 0;
            d2N_xi_xi_k_I(1, 1, :, 18) = 0;
            d2N_xi_xi_k_I(1, 1, :, 19) = 0;
            d2N_xi_xi_k_I(1, 1, :, 20) = 0;
            d2N_xi_xi_k_I(1, 2, :, 1) = (1 / 8) * (zeta .* zeta + 2 * xi .* zeta + 2 * eta .* zeta - 2 * xi - 2 * eta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 2) = ((xi + 1) .* (zeta - 1)) / 8 - ((zeta - 1) .* (eta - xi + zeta + 2)) / 8 - ((eta - 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 3) = -((xi + 1) .* (zeta - 1)) / 8 - ((zeta - 1) .* (eta + xi - zeta - 2)) / 8 - ((eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 4) = ((eta + 1) .* (zeta - 1)) / 8 - ((zeta - 1) .* (xi - eta + zeta + 2)) / 8 - ((xi - 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 5) = (1 / 8) * (zeta .* zeta - 2 * xi .* zeta - 2 * eta .* zeta - 2 * xi - 2 * eta + zeta);
            d2N_xi_xi_k_I(1, 2, :, 6) = (1 / 8) * (-zeta .* zeta - 2 * xi .* zeta + 2 * eta .* zeta - 2 * xi + 2 * eta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 7) = (1 / 8) * (zeta .* zeta + 2 * xi .* zeta + 2 * eta .* zeta + 2 * xi + 2 * eta + zeta);
            d2N_xi_xi_k_I(1, 2, :, 8) = (1 / 8) * (-zeta .* zeta + 2 * xi .* zeta - 2 * eta .* zeta + 2 * xi - 2 * eta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 9) = 0.5 * (-xi .* zeta + xi);
            d2N_xi_xi_k_I(1, 2, :, 10) = 0.5 * (eta .* zeta - eta);
            d2N_xi_xi_k_I(1, 2, :, 11) = 0.5 * (xi .* zeta - xi);
            d2N_xi_xi_k_I(1, 2, :, 12) = 0.5 * (-eta .* zeta + eta);
            d2N_xi_xi_k_I(1, 2, :, 13) = (xi .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(1, 2, :, 14) = -(eta .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(1, 2, :, 15) = 0.5 * (-xi .* zeta - xi);
            d2N_xi_xi_k_I(1, 2, :, 16) = 0.5 * (eta .* zeta + eta);
            d2N_xi_xi_k_I(1, 2, :, 17) = 0.25 * (-zeta .* zeta + 1);
            d2N_xi_xi_k_I(1, 2, :, 18) = 0.25 * (zeta .* zeta - 1);
            d2N_xi_xi_k_I(1, 2, :, 19) = 1 / 4 - zeta.^2 / 4;
            d2N_xi_xi_k_I(1, 2, :, 20) = 0.25 * (zeta .* zeta - 1);
            d2N_xi_xi_k_I(1, 3, :, 1) = (1 / 8) * (eta .* eta + 2 * xi .* eta + 2 * eta .* zeta - 2 * xi - eta - 2 * zeta);
            d2N_xi_xi_k_I(1, 3, :, 2) = ((eta - 1) .* (xi + 1)) / 8 - ((eta - 1) .* (eta - xi + zeta + 2)) / 8 - ((eta - 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 3, :, 3) = ((eta + 1) .* (zeta - 1)) / 8 - ((eta + 1) .* (xi + 1)) / 8 - ((eta + 1) .* (eta + xi - zeta - 2)) / 8;
            d2N_xi_xi_k_I(1, 3, :, 4) = -((eta + 1) .* (xi - eta + zeta + 2)) / 8 - ((eta + 1) .* (xi - 1)) / 8 - ((eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 3, :, 5) = (1 / 8) * (-eta .* eta - 2 * xi .* eta + 2 * eta .* zeta + 2 * xi + eta - 2 * zeta);
            d2N_xi_xi_k_I(1, 3, :, 6) = (1 / 8) * (eta .* eta - 2 * xi .* eta - 2 * eta .* zeta + 2 * xi - eta + 2 * zeta);
            d2N_xi_xi_k_I(1, 3, :, 7) = (1 / 8) * (eta .* eta + 2 * xi .* eta + 2 * eta .* zeta + 2 * xi + eta + 2 * zeta);
            d2N_xi_xi_k_I(1, 3, :, 8) = (1 / 8) * (-eta .* eta + 2 * xi .* eta - 2 * eta .* zeta + 2 * xi - eta - 2 * zeta);
            d2N_xi_xi_k_I(1, 3, :, 9) = 0.5 * (-xi .* eta + xi);
            d2N_xi_xi_k_I(1, 3, :, 10) = 0.25 * (eta .* eta - 1);
            d2N_xi_xi_k_I(1, 3, :, 11) = 0.5 * (xi .* eta + xi);
            d2N_xi_xi_k_I(1, 3, :, 12) = 0.25 * (-eta .* eta + 1);
            d2N_xi_xi_k_I(1, 3, :, 13) = (xi .* (eta - 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 14) = 1 / 4 - eta.^2 / 4;
            d2N_xi_xi_k_I(1, 3, :, 15) = 0.5 * (-xi .* eta - xi);
            d2N_xi_xi_k_I(1, 3, :, 16) = 0.25 * (eta .* eta - 1);
            d2N_xi_xi_k_I(1, 3, :, 17) = 0.5 * (-eta .* zeta + zeta);
            d2N_xi_xi_k_I(1, 3, :, 18) = 0.5 * (eta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 3, :, 19) = 0.5 * (-eta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 3, :, 20) = 0.5 * (eta .* zeta + zeta);
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = d2N_xi_xi_k_I(1, 2, :, 1);
            d2N_xi_xi_k_I(2, 1, :, 2) = d2N_xi_xi_k_I(1, 2, :, 2);
            d2N_xi_xi_k_I(2, 1, :, 3) = d2N_xi_xi_k_I(1, 2, :, 3);
            d2N_xi_xi_k_I(2, 1, :, 4) = d2N_xi_xi_k_I(1, 2, :, 4);
            d2N_xi_xi_k_I(2, 1, :, 5) = d2N_xi_xi_k_I(1, 2, :, 5);
            d2N_xi_xi_k_I(2, 1, :, 6) = d2N_xi_xi_k_I(1, 2, :, 6);
            d2N_xi_xi_k_I(2, 1, :, 7) = d2N_xi_xi_k_I(1, 2, :, 7);
            d2N_xi_xi_k_I(2, 1, :, 8) = d2N_xi_xi_k_I(1, 2, :, 8);
            d2N_xi_xi_k_I(2, 1, :, 9) = d2N_xi_xi_k_I(1, 2, :, 9);
            d2N_xi_xi_k_I(2, 1, :, 10) = d2N_xi_xi_k_I(1, 2, :, 10);
            d2N_xi_xi_k_I(2, 1, :, 11) = d2N_xi_xi_k_I(1, 2, :, 11);
            d2N_xi_xi_k_I(2, 1, :, 12) = d2N_xi_xi_k_I(1, 2, :, 12);
            d2N_xi_xi_k_I(2, 1, :, 13) = d2N_xi_xi_k_I(1, 2, :, 13);
            d2N_xi_xi_k_I(2, 1, :, 14) = d2N_xi_xi_k_I(1, 2, :, 14);
            d2N_xi_xi_k_I(2, 1, :, 15) = d2N_xi_xi_k_I(1, 2, :, 15);
            d2N_xi_xi_k_I(2, 1, :, 16) = d2N_xi_xi_k_I(1, 2, :, 16);
            d2N_xi_xi_k_I(2, 1, :, 17) = d2N_xi_xi_k_I(1, 2, :, 17);
            d2N_xi_xi_k_I(2, 1, :, 18) = d2N_xi_xi_k_I(1, 2, :, 18);
            d2N_xi_xi_k_I(2, 1, :, 19) = d2N_xi_xi_k_I(1, 2, :, 19);
            d2N_xi_xi_k_I(2, 1, :, 20) = d2N_xi_xi_k_I(1, 2, :, 20);
            d2N_xi_xi_k_I(2, 2, :, 1) = 0.25 * (xi .* zeta - xi - zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 2) = 0.25 * (-xi .* zeta + xi - zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 3) = 0.25 * (-xi .* zeta + xi - zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 4) = 0.25 * (xi .* zeta - xi - zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 5) = 0.25 * (-xi .* zeta - xi + zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 6) = 0.25 * (xi .* zeta + xi + zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 7) = 0.25 * (xi .* zeta + xi + zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 8) = 0.25 * (-xi .* zeta - xi + zeta + 1);
            d2N_xi_xi_k_I(2, 2, :, 9) = 0;
            d2N_xi_xi_k_I(2, 2, :, 10) = ((xi + 1) .* (zeta - 1)) / 2;
            d2N_xi_xi_k_I(2, 2, :, 11) = 0;
            d2N_xi_xi_k_I(2, 2, :, 12) = 0.5 * (-xi .* zeta + xi + zeta - 1);
            d2N_xi_xi_k_I(2, 2, :, 13) = 0;
            d2N_xi_xi_k_I(2, 2, :, 14) = 0.5 * (-xi .* zeta - xi - zeta - 1);
            d2N_xi_xi_k_I(2, 2, :, 15) = 0;
            d2N_xi_xi_k_I(2, 2, :, 16) = 0.5 * (xi .* zeta + xi - zeta - 1);
            d2N_xi_xi_k_I(2, 2, :, 17) = 0;
            d2N_xi_xi_k_I(2, 2, :, 18) = 0;
            d2N_xi_xi_k_I(2, 2, :, 19) = 0;
            d2N_xi_xi_k_I(2, 2, :, 20) = 0;
            d2N_xi_xi_k_I(2, 3, :, 1) = ((xi - 1) .* (zeta - 1)) / 8 + ((xi - 1) .* (eta + xi + zeta + 2)) / 8 + ((eta - 1) .* (xi - 1)) / 8;
            d2N_xi_xi_k_I(2, 3, :, 2) = -((xi + 1) .* (zeta - 1)) / 8 - ((xi + 1) .* (eta - xi + zeta + 2)) / 8 - ((eta - 1) .* (xi + 1)) / 8;
            d2N_xi_xi_k_I(2, 3, :, 3) = ((xi + 1) .* (zeta - 1)) / 8 - ((xi + 1) .* (eta + xi - zeta - 2)) / 8 - ((eta + 1) .* (xi + 1)) / 8;
            d2N_xi_xi_k_I(2, 3, :, 4) = (1 / 8) * (-xi .* xi + 2 * xi .* eta - 2 * xi .* zeta + xi - 2 * eta + 2 * zeta);
            d2N_xi_xi_k_I(2, 3, :, 5) = (1 / 8) * (-xi .* xi - 2 * xi .* eta + 2 * xi .* zeta + xi + 2 * eta - 2 * zeta);
            d2N_xi_xi_k_I(2, 3, :, 6) = (1 / 8) * (-xi .* xi + 2 * xi .* eta - 2 * xi .* zeta - xi + 2 * eta - 2 * zeta);
            d2N_xi_xi_k_I(2, 3, :, 7) = ((xi + 1) .* (zeta + 1)) / 8 + ((xi + 1) .* (eta + xi + zeta - 2)) / 8 + ((eta + 1) .* (xi + 1)) / 8;
            d2N_xi_xi_k_I(2, 3, :, 8) = (1 / 8) * (xi .* xi - 2 * xi .* eta - 2 * xi .* zeta - xi + 2 * eta + 2 * zeta);
            d2N_xi_xi_k_I(2, 3, :, 9) = 0.25 * (-xi .* xi + 1);
            d2N_xi_xi_k_I(2, 3, :, 10) = (eta .* (xi + 1)) / 2;
            d2N_xi_xi_k_I(2, 3, :, 11) = xi.^2 / 4 - 1 / 4;
            d2N_xi_xi_k_I(2, 3, :, 12) = 0.5 * (-xi .* eta + eta);
            d2N_xi_xi_k_I(2, 3, :, 13) = 0.25 * (xi .* xi - 1);
            d2N_xi_xi_k_I(2, 3, :, 14) = 0.5 * (-xi .* eta - eta);
            d2N_xi_xi_k_I(2, 3, :, 15) = 0.25 * (-xi .* xi + 1);
            d2N_xi_xi_k_I(2, 3, :, 16) = 0.5 * (xi .* eta - eta);
            d2N_xi_xi_k_I(2, 3, :, 17) = -(zeta .* (xi - 1)) / 2;
            d2N_xi_xi_k_I(2, 3, :, 18) = (zeta .* (xi + 1)) / 2;
            d2N_xi_xi_k_I(2, 3, :, 19) = -(zeta .* (xi + 1)) / 2;
            d2N_xi_xi_k_I(2, 3, :, 20) = (zeta .* (xi - 1)) / 2;
            % second derivative wrt zeta
            d2N_xi_xi_k_I(3, 1, :, 1) = d2N_xi_xi_k_I(1, 3, :, 1);
            d2N_xi_xi_k_I(3, 1, :, 2) = d2N_xi_xi_k_I(1, 3, :, 2);
            d2N_xi_xi_k_I(3, 1, :, 3) = d2N_xi_xi_k_I(1, 3, :, 3);
            d2N_xi_xi_k_I(3, 1, :, 4) = d2N_xi_xi_k_I(1, 3, :, 4);
            d2N_xi_xi_k_I(3, 1, :, 5) = d2N_xi_xi_k_I(1, 3, :, 5);
            d2N_xi_xi_k_I(3, 1, :, 6) = d2N_xi_xi_k_I(1, 3, :, 6);
            d2N_xi_xi_k_I(3, 1, :, 7) = d2N_xi_xi_k_I(1, 3, :, 7);
            d2N_xi_xi_k_I(3, 1, :, 8) = d2N_xi_xi_k_I(1, 3, :, 8);
            d2N_xi_xi_k_I(3, 1, :, 9) = d2N_xi_xi_k_I(1, 3, :, 9);
            d2N_xi_xi_k_I(3, 1, :, 10) = d2N_xi_xi_k_I(1, 3, :, 10);
            d2N_xi_xi_k_I(3, 1, :, 11) = d2N_xi_xi_k_I(1, 3, :, 11);
            d2N_xi_xi_k_I(3, 1, :, 12) = d2N_xi_xi_k_I(1, 3, :, 12);
            d2N_xi_xi_k_I(3, 1, :, 13) = d2N_xi_xi_k_I(1, 3, :, 13);
            d2N_xi_xi_k_I(3, 1, :, 14) = d2N_xi_xi_k_I(1, 3, :, 14);
            d2N_xi_xi_k_I(3, 1, :, 15) = d2N_xi_xi_k_I(1, 3, :, 15);
            d2N_xi_xi_k_I(3, 1, :, 16) = d2N_xi_xi_k_I(1, 3, :, 16);
            d2N_xi_xi_k_I(3, 1, :, 17) = d2N_xi_xi_k_I(1, 3, :, 17);
            d2N_xi_xi_k_I(3, 1, :, 18) = d2N_xi_xi_k_I(1, 3, :, 18);
            d2N_xi_xi_k_I(3, 1, :, 19) = d2N_xi_xi_k_I(1, 3, :, 19);
            d2N_xi_xi_k_I(3, 1, :, 20) = d2N_xi_xi_k_I(1, 3, :, 20);
            d2N_xi_xi_k_I(3, 2, :, 1) = d2N_xi_xi_k_I(2, 3, :, 1);
            d2N_xi_xi_k_I(3, 2, :, 2) = d2N_xi_xi_k_I(2, 3, :, 2);
            d2N_xi_xi_k_I(3, 2, :, 3) = d2N_xi_xi_k_I(2, 3, :, 3);
            d2N_xi_xi_k_I(3, 2, :, 4) = d2N_xi_xi_k_I(2, 3, :, 4);
            d2N_xi_xi_k_I(3, 2, :, 5) = d2N_xi_xi_k_I(2, 3, :, 5);
            d2N_xi_xi_k_I(3, 2, :, 6) = d2N_xi_xi_k_I(2, 3, :, 6);
            d2N_xi_xi_k_I(3, 2, :, 7) = d2N_xi_xi_k_I(2, 3, :, 7);
            d2N_xi_xi_k_I(3, 2, :, 8) = d2N_xi_xi_k_I(2, 3, :, 8);
            d2N_xi_xi_k_I(3, 2, :, 9) = d2N_xi_xi_k_I(2, 3, :, 9);
            d2N_xi_xi_k_I(3, 2, :, 10) = d2N_xi_xi_k_I(2, 3, :, 10);
            d2N_xi_xi_k_I(3, 2, :, 11) = d2N_xi_xi_k_I(2, 3, :, 11);
            d2N_xi_xi_k_I(3, 2, :, 12) = d2N_xi_xi_k_I(2, 3, :, 12);
            d2N_xi_xi_k_I(3, 2, :, 13) = d2N_xi_xi_k_I(2, 3, :, 13);
            d2N_xi_xi_k_I(3, 2, :, 14) = d2N_xi_xi_k_I(2, 3, :, 14);
            d2N_xi_xi_k_I(3, 2, :, 15) = d2N_xi_xi_k_I(2, 3, :, 15);
            d2N_xi_xi_k_I(3, 2, :, 16) = d2N_xi_xi_k_I(2, 3, :, 16);
            d2N_xi_xi_k_I(3, 2, :, 17) = d2N_xi_xi_k_I(2, 3, :, 17);
            d2N_xi_xi_k_I(3, 2, :, 18) = d2N_xi_xi_k_I(2, 3, :, 18);
            d2N_xi_xi_k_I(3, 2, :, 19) = d2N_xi_xi_k_I(2, 3, :, 19);
            d2N_xi_xi_k_I(3, 2, :, 20) = d2N_xi_xi_k_I(2, 3, :, 20);
            d2N_xi_xi_k_I(3, 3, :, 1) = 0.25 * (xi .* eta - xi - eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 2) = 0.25 * (-xi .* eta + xi - eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 3) = 0.25 * (xi .* eta + xi + eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 4) = 0.25 * (-xi .* eta - xi + eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 5) = 0.25 * (xi .* eta - xi - eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 6) = 0.25 * (-xi .* eta + xi - eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 7) = 0.25 * (xi .* eta + xi + eta + 1);
            d2N_xi_xi_k_I(3, 3, :, 8) = -((eta + 1) .* (xi - 1)) / 4;
            d2N_xi_xi_k_I(3, 3, :, 9) = 0;
            d2N_xi_xi_k_I(3, 3, :, 10) = 0;
            d2N_xi_xi_k_I(3, 3, :, 11) = 0;
            d2N_xi_xi_k_I(3, 3, :, 12) = 0;
            d2N_xi_xi_k_I(3, 3, :, 13) = 0;
            d2N_xi_xi_k_I(3, 3, :, 14) = 0;
            d2N_xi_xi_k_I(3, 3, :, 15) = 0;
            d2N_xi_xi_k_I(3, 3, :, 16) = 0;
            d2N_xi_xi_k_I(3, 3, :, 17) = 0.5 * (-xi .* eta + xi + eta - 1);
            d2N_xi_xi_k_I(3, 3, :, 18) = 0.5 * (xi .* eta - xi + eta - 1);
            d2N_xi_xi_k_I(3, 3, :, 19) = 0.5 * (-xi .* eta - xi - eta - 1);
            d2N_xi_xi_k_I(3, 3, :, 20) = 0.5 * (xi .* eta + xi - eta - 1);
        elseif numberOfNodes == 27
            N_k_I(:, 1) = xi .* eta .* zeta .* (xi - 1) .* (eta - 1) .* (zeta - 1) / 8;
            N_k_I(:, 2) = xi .* eta .* zeta .* (xi + 1) .* (eta - 1) .* (zeta - 1) / 8;
            N_k_I(:, 3) = xi .* eta .* zeta .* (xi + 1) .* (eta + 1) .* (zeta - 1) / 8;
            N_k_I(:, 4) = xi .* eta .* zeta .* (xi - 1) .* (eta + 1) .* (zeta - 1) / 8;
            N_k_I(:, 5) = xi .* eta .* zeta .* (xi - 1) .* (eta - 1) .* (zeta + 1) / 8;
            N_k_I(:, 6) = xi .* eta .* zeta .* (xi + 1) .* (eta - 1) .* (zeta + 1) / 8;
            N_k_I(:, 7) = xi .* eta .* zeta .* (xi + 1) .* (eta + 1) .* (zeta + 1) / 8;
            N_k_I(:, 8) = xi .* eta .* zeta .* (xi - 1) .* (eta + 1) .* (zeta + 1) / 8;
            N_k_I(:, 9) = xi .* eta .* (xi - 1) .* (eta - 1) .* (1 - zeta.^2) / 4;
            N_k_I(:, 10) = xi .* eta .* (xi + 1) .* (eta - 1) .* (1 - zeta.^2) / 4;
            N_k_I(:, 11) = xi .* eta .* (xi + 1) .* (eta + 1) .* (1 - zeta.^2) / 4;
            N_k_I(:, 12) = xi .* eta .* (xi - 1) .* (eta + 1) .* (1 - zeta.^2) / 4;
            N_k_I(:, 13) = eta .* zeta .* (1 - xi.^2) .* (eta - 1) .* (zeta - 1) / 4;
            N_k_I(:, 14) = xi .* zeta .* (xi + 1) .* (1 - eta.^2) .* (zeta - 1) / 4;
            N_k_I(:, 15) = eta .* zeta .* (1 - xi.^2) .* (eta + 1) .* (zeta - 1) / 4;
            N_k_I(:, 16) = xi .* zeta .* (xi - 1) .* (1 - eta.^2) .* (zeta - 1) / 4;
            N_k_I(:, 17) = zeta .* (1 - xi.^2) .* (1 - eta.^2) .* (zeta - 1) / 2;
            N_k_I(:, 18) = eta .* zeta .* (1 - xi.^2) .* (eta - 1) .* (zeta + 1) / 4;
            N_k_I(:, 19) = xi .* zeta .* (xi + 1) .* (1 - eta.^2) .* (zeta + 1) / 4;
            N_k_I(:, 20) = eta .* zeta .* (1 - xi.^2) .* (eta + 1) .* (zeta + 1) / 4;
            N_k_I(:, 21) = xi .* zeta .* (xi - 1) .* (1 - eta.^2) .* (zeta + 1) / 4;
            N_k_I(:, 22) = zeta .* (1 - xi.^2) .* (1 - eta.^2) .* (zeta + 1) / 2;
            N_k_I(:, 23) = eta .* (1 - xi.^2) .* (eta - 1) .* (1 - zeta.^2) / 2;
            N_k_I(:, 24) = xi .* (xi + 1) .* (1 - eta.^2) .* (1 - zeta.^2) / 2;
            N_k_I(:, 25) = eta .* (1 - xi.^2) .* (eta + 1) .* (1 - zeta.^2) / 2;
            N_k_I(:, 26) = xi .* (xi - 1) .* (1 - eta.^2) .* (1 - zeta.^2) / 2;
            N_k_I(:, 27) = (1 - xi.^2) .* (1 - eta.^2) .* (1 - zeta.^2);
            % derivative wrt xi
            dN_xi_k_I(1, :, 1) = (eta .* zeta .* (eta - 1) .* (zeta - 1) .* (2 * xi - 1)) / 8;
            dN_xi_k_I(1, :, 2) = (eta .* zeta .* (eta - 1) .* (zeta - 1) .* (2 * xi + 1)) / 8;
            dN_xi_k_I(1, :, 3) = (eta .* zeta .* (eta + 1) .* (zeta - 1) .* (2 * xi + 1)) / 8;
            dN_xi_k_I(1, :, 4) = (eta .* zeta .* (eta + 1) .* (zeta - 1) .* (2 * xi - 1)) / 8;
            dN_xi_k_I(1, :, 5) = (eta .* zeta .* (eta - 1) .* (zeta + 1) .* (2 * xi - 1)) / 8;
            dN_xi_k_I(1, :, 6) = (eta .* zeta .* (eta - 1) .* (zeta + 1) .* (2 * xi + 1)) / 8;
            dN_xi_k_I(1, :, 7) = (eta .* zeta .* (eta + 1) .* (zeta + 1) .* (2 * xi + 1)) / 8;
            dN_xi_k_I(1, :, 8) = (eta .* zeta .* (eta + 1) .* (zeta + 1) .* (2 * xi - 1)) / 8;
            dN_xi_k_I(1, :, 9) = (eta .* (eta - 1) .* (1 - zeta.^2) .* (2 * xi - 1)) / 4;
            dN_xi_k_I(1, :, 10) = (eta .* (eta - 1) .* (1 - zeta.^2) .* (2 * xi + 1)) / 4;
            dN_xi_k_I(1, :, 11) = (eta .* (eta + 1) .* (1 - zeta.^2) .* (2 * xi + 1)) / 4;
            dN_xi_k_I(1, :, 12) = (eta .* (eta + 1) .* (1 - zeta.^2) .* (2 * xi - 1)) / 4;
            dN_xi_k_I(1, :, 13) = (eta .* zeta .* (eta - 1) .* (zeta - 1) .* (-2 * xi)) / 4;
            dN_xi_k_I(1, :, 14) = (zeta .* (1 - eta.^2) .* (zeta - 1) .* (2 * xi + 1)) / 4;
            dN_xi_k_I(1, :, 15) = (eta .* zeta .* (eta + 1) .* (zeta - 1) .* (-2 * xi)) / 4;
            dN_xi_k_I(1, :, 16) = (zeta .* (1 - eta.^2) .* (zeta - 1) .* (2 * xi - 1)) / 4;
            dN_xi_k_I(1, :, 17) = (zeta .* (1 - eta.^2) .* (zeta - 1) .* (-2 * xi)) / 2;
            dN_xi_k_I(1, :, 18) = (eta .* zeta .* (eta - 1) .* (zeta + 1) .* (-2 * xi)) / 4;
            dN_xi_k_I(1, :, 19) = (zeta .* (1 - eta.^2) .* (zeta + 1) .* (2 * xi + 1)) / 4;
            dN_xi_k_I(1, :, 20) = (eta .* zeta .* (eta + 1) .* (zeta + 1) .* (-2 * xi)) / 4;
            dN_xi_k_I(1, :, 21) = (zeta .* (1 - eta.^2) .* (zeta + 1) .* (2 * xi - 1)) / 4;
            dN_xi_k_I(1, :, 22) = (zeta .* (1 - eta.^2) .* (zeta + 1) .* (-2 * xi)) / 2;
            dN_xi_k_I(1, :, 23) = (eta .* (eta - 1) .* (1 - zeta.^2) .* (-2 * xi)) / 2;
            dN_xi_k_I(1, :, 24) = ((1 - eta.^2) .* (1 - zeta.^2) .* (2 * xi + 1)) / 2;
            dN_xi_k_I(1, :, 25) = (eta .* (eta + 1) .* (1 - zeta.^2) .* (-2 * xi)) / 2;
            dN_xi_k_I(1, :, 26) = ((1 - eta.^2) .* (1 - zeta.^2) .* (2 * xi - 1)) / 2;
            dN_xi_k_I(1, :, 27) = ((1 - eta.^2) .* (1 - zeta.^2) .* (-2 * xi));
            % derivative wrt eta
            dN_xi_k_I(2, :, 1) = (xi .* zeta .* (xi - 1) .* (zeta - 1) .* (2 * eta - 1)) / 8;
            dN_xi_k_I(2, :, 2) = (xi .* zeta .* (xi + 1) .* (zeta - 1) .* (2 * eta - 1)) / 8;
            dN_xi_k_I(2, :, 3) = (xi .* zeta .* (xi + 1) .* (zeta - 1) .* (2 * eta + 1)) / 8;
            dN_xi_k_I(2, :, 4) = (xi .* zeta .* (xi - 1) .* (zeta - 1) .* (2 * eta + 1)) / 8;
            dN_xi_k_I(2, :, 5) = (xi .* zeta .* (xi - 1) .* (zeta + 1) .* (2 * eta - 1)) / 8;
            dN_xi_k_I(2, :, 6) = (xi .* zeta .* (xi + 1) .* (zeta + 1) .* (2 * eta - 1)) / 8;
            dN_xi_k_I(2, :, 7) = (xi .* zeta .* (xi + 1) .* (zeta + 1) .* (2 * eta + 1)) / 8;
            dN_xi_k_I(2, :, 8) = (xi .* zeta .* (xi - 1) .* (zeta + 1) .* (2 * eta + 1)) / 8;
            dN_xi_k_I(2, :, 9) = (xi .* (xi - 1) .* (1 - zeta.^2) .* (2 * eta - 1)) / 4;
            dN_xi_k_I(2, :, 10) = (xi .* (xi + 1) .* (1 - zeta.^2) .* (2 * eta - 1)) / 4;
            dN_xi_k_I(2, :, 11) = (xi .* (xi + 1) .* (1 - zeta.^2) .* (2 * eta + 1)) / 4;
            dN_xi_k_I(2, :, 12) = (xi .* (xi - 1) .* (1 - zeta.^2) .* (2 * eta + 1)) / 4;
            dN_xi_k_I(2, :, 13) = (zeta .* (1 - xi.^2) .* (zeta - 1) .* (2 * eta - 1)) / 4;
            dN_xi_k_I(2, :, 14) = (xi .* zeta .* (xi + 1) .* (zeta - 1) .* (-2 * eta)) / 4;
            dN_xi_k_I(2, :, 15) = (zeta .* (1 - xi.^2) .* (zeta - 1) .* (2 * eta + 1)) / 4;
            dN_xi_k_I(2, :, 16) = (xi .* zeta .* (xi - 1) .* (zeta - 1) .* (-2 * eta)) / 4;
            dN_xi_k_I(2, :, 17) = (zeta .* (1 - xi.^2) .* (zeta - 1) .* (-2 * eta)) / 2;
            dN_xi_k_I(2, :, 18) = (zeta .* (1 - xi.^2) .* (zeta + 1) .* (2 * eta - 1)) / 4;
            dN_xi_k_I(2, :, 19) = (xi .* zeta .* (xi + 1) .* (zeta + 1) .* (-2 * eta)) / 4;
            dN_xi_k_I(2, :, 20) = (zeta .* (1 - xi.^2) .* (zeta + 1) .* (2 * eta + 1)) / 4;
            dN_xi_k_I(2, :, 21) = (xi .* zeta .* (xi - 1) .* (zeta + 1) .* (-2 * eta)) / 4;
            dN_xi_k_I(2, :, 22) = (zeta .* (1 - xi.^2) .* (zeta + 1) .* (-2 * eta)) / 2;
            dN_xi_k_I(2, :, 23) = ((1 - xi.^2) .* (1 - zeta.^2) .* (2 * eta - 1)) / 2;
            dN_xi_k_I(2, :, 24) = (xi .* (xi + 1) .* (1 - zeta.^2) .* (-2 * eta)) / 2;
            dN_xi_k_I(2, :, 25) = ((1 - xi.^2) .* (1 - zeta.^2) .* (2 * eta + 1)) / 2;
            dN_xi_k_I(2, :, 26) = (xi .* (xi - 1) .* (1 - zeta.^2) .* (-2 * eta)) / 2;
            dN_xi_k_I(2, :, 27) = ((1 - xi.^2) .* (1 - zeta.^2) .* (-2 * eta));
            % derivative wrt zeta
            dN_xi_k_I(3, :, 1) = (xi .* eta .* (xi - 1) .* (eta - 1) .* (2 * zeta - 1)) / 8;
            dN_xi_k_I(3, :, 2) = (xi .* eta .* (xi + 1) .* (eta - 1) .* (2 * zeta - 1)) / 8;
            dN_xi_k_I(3, :, 3) = (xi .* eta .* (xi + 1) .* (eta + 1) .* (2 * zeta - 1)) / 8;
            dN_xi_k_I(3, :, 4) = (xi .* eta .* (xi - 1) .* (eta + 1) .* (2 * zeta - 1)) / 8;
            dN_xi_k_I(3, :, 5) = (xi .* eta .* (xi - 1) .* (eta - 1) .* (2 * zeta + 1)) / 8;
            dN_xi_k_I(3, :, 6) = (xi .* eta .* (xi + 1) .* (eta - 1) .* (2 * zeta + 1)) / 8;
            dN_xi_k_I(3, :, 7) = (xi .* eta .* (xi + 1) .* (eta + 1) .* (2 * zeta + 1)) / 8;
            dN_xi_k_I(3, :, 8) = (xi .* eta .* (xi - 1) .* (eta + 1) .* (2 * zeta + 1)) / 8;
            dN_xi_k_I(3, :, 9) = (xi .* eta .* (xi - 1) .* (eta - 1) .* (-2 * zeta)) / 4;
            dN_xi_k_I(3, :, 10) = (xi .* eta .* (xi + 1) .* (eta - 1) .* (-2 * zeta)) / 4;
            dN_xi_k_I(3, :, 11) = (xi .* eta .* (xi + 1) .* (eta + 1) .* (-2 * zeta)) / 4;
            dN_xi_k_I(3, :, 12) = (xi .* eta .* (xi - 1) .* (eta + 1) .* (-2 * zeta)) / 4;
            dN_xi_k_I(3, :, 13) = (eta .* (1 - xi.^2) .* (eta - 1) .* (2 * zeta - 1)) / 4;
            dN_xi_k_I(3, :, 14) = (xi .* (xi + 1) .* (1 - eta.^2) .* (2 * zeta - 1)) / 4;
            dN_xi_k_I(3, :, 15) = (eta .* (1 - xi.^2) .* (eta + 1) .* (2 * zeta - 1)) / 4;
            dN_xi_k_I(3, :, 16) = (xi .* (xi - 1) .* (1 - eta.^2) .* (2 * zeta - 1)) / 4;
            dN_xi_k_I(3, :, 17) = ((1 - xi.^2) .* (1 - eta.^2) .* (2 * zeta - 1)) / 2;
            dN_xi_k_I(3, :, 18) = (eta .* (1 - xi.^2) .* (eta - 1) .* (2 * zeta + 1)) / 4;
            dN_xi_k_I(3, :, 19) = (xi .* (xi + 1) .* (1 - eta.^2) .* (2 * zeta + 1)) / 4;
            dN_xi_k_I(3, :, 20) = (eta .* (1 - xi.^2) .* (eta + 1) .* (2 * zeta + 1)) / 4;
            dN_xi_k_I(3, :, 21) = (xi .* (xi - 1) .* (1 - eta.^2) .* (2 * zeta + 1)) / 4;
            dN_xi_k_I(3, :, 22) = ((1 - xi.^2) .* (1 - eta.^2) .* (2 * zeta + 1)) / 2;
            dN_xi_k_I(3, :, 23) = (eta .* (1 - xi.^2) .* (eta - 1) .* (-2 * zeta)) / 2;
            dN_xi_k_I(3, :, 24) = (xi .* (xi + 1) .* (1 - eta.^2) .* (-2 * zeta)) / 2;
            dN_xi_k_I(3, :, 25) = (eta .* (1 - xi.^2) .* (eta + 1) .* (-2 * zeta)) / 2;
            dN_xi_k_I(3, :, 26) = (xi .* (xi - 1) .* (1 - eta.^2) .* (-2 * zeta)) / 2;
            dN_xi_k_I(3, :, 27) = ((1 - xi.^2) .* (1 - eta.^2) .* (-2 * zeta));
            % second derivative wrt xi
            d2N_xi_xi_k_I(1, 1, :, 1) = (1 / 4) * (eta .* eta .* zeta .* zeta - eta .* zeta .* zeta - eta .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 2) = (1 / 4) * (eta .* eta .* zeta .* zeta - eta .* zeta .* zeta - eta .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 3) = (1 / 4) * (eta .* eta .* zeta .* zeta + eta .* zeta .* zeta - eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 4) = (1 / 4) * (eta .* eta .* zeta .* zeta + eta .* zeta .* zeta - eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 5) = (1 / 4) * (eta .* eta .* zeta .* zeta - eta .* zeta .* zeta + eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 6) = (1 / 4) * (eta .* eta .* zeta .* zeta - eta .* zeta .* zeta + eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 7) = (eta .* zeta .* (eta + 1) .* (zeta + 1)) / 4;
            d2N_xi_xi_k_I(1, 1, :, 8) = (1 / 4) * (eta .* eta .* zeta .* zeta + eta .* zeta .* zeta + eta .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 9) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* zeta .* zeta + eta .* eta - eta);
            d2N_xi_xi_k_I(1, 1, :, 10) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* zeta .* zeta + eta .* eta - eta);
            d2N_xi_xi_k_I(1, 1, :, 11) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* zeta .* zeta + eta .* eta + eta);
            d2N_xi_xi_k_I(1, 1, :, 12) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* zeta .* zeta + eta .* eta + eta);
            d2N_xi_xi_k_I(1, 1, :, 13) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* zeta .* zeta + eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 14) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* eta .* zeta + zeta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 1, :, 15) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* zeta .* zeta + eta .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 16) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* eta .* zeta + zeta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 1, :, 17) = eta .* eta .* zeta .* zeta - eta .* eta .* zeta - zeta .* zeta + zeta;
            d2N_xi_xi_k_I(1, 1, :, 18) = 0.5 * (-eta .* eta .* zeta .* zeta + eta .* zeta .* zeta - eta .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 19) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* eta .* zeta + zeta .* zeta + zeta);
            d2N_xi_xi_k_I(1, 1, :, 20) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* zeta .* zeta - eta .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 1, :, 21) = 0.5 * (-eta .* eta .* zeta .* zeta - eta .* eta .* zeta + zeta .* zeta + zeta);
            d2N_xi_xi_k_I(1, 1, :, 22) = eta .* eta .* zeta .* zeta + eta .* eta .* zeta - zeta .* zeta - zeta;
            d2N_xi_xi_k_I(1, 1, :, 23) = eta .* eta .* zeta .* zeta - eta .* zeta .* zeta - eta .* eta + eta;
            d2N_xi_xi_k_I(1, 1, :, 24) = eta .* eta .* zeta .* zeta - eta .* eta - zeta .* zeta + 1;
            d2N_xi_xi_k_I(1, 1, :, 25) = eta .* eta .* zeta .* zeta + eta .* zeta .* zeta - eta .* eta - eta;
            d2N_xi_xi_k_I(1, 1, :, 26) = eta .* eta .* zeta .* zeta - eta .* eta - zeta .* zeta + 1;
            d2N_xi_xi_k_I(1, 1, :, 27) = -2 * eta .* eta .* zeta .* zeta + 2 * eta .* eta + 2 * zeta .* zeta - 2;
            d2N_xi_xi_k_I(1, 2, :, 1) = (1 / 8) * (4 * xi .* eta .* zeta .* zeta - 2 * xi .* zeta .* zeta - 2 .* eta .* zeta .* zeta - 4 * xi .* eta .* zeta + zeta .* zeta + 2 * xi .* zeta + 2 * eta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 2) = (1 / 8) * (4 * xi .* eta .* zeta .* zeta - 2 * xi .* zeta .* zeta + 2 .* eta .* zeta .* zeta - 4 * xi .* eta .* zeta - zeta .* zeta + 2 * xi .* zeta - 2 * eta .* zeta + zeta);
            d2N_xi_xi_k_I(1, 2, :, 3) = (eta .* zeta .* (2 * xi + 1) .* (zeta - 1)) / 8 + (zeta .* (2 * xi + 1) .* (eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 4) = (eta .* zeta .* (2 * xi - 1) .* (zeta - 1)) / 8 + (zeta .* (2 * xi - 1) .* (eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 5) = (1 / 8) * (4 * xi .* eta .* zeta .* zeta - 2 * xi .* zeta .* zeta - 2 .* eta .* zeta .* zeta + 4 * xi .* eta .* zeta + zeta .* zeta - 2 * xi .* zeta - 2 * eta .* zeta + zeta);
            d2N_xi_xi_k_I(1, 2, :, 6) = (1 / 8) * (4 * xi .* eta .* zeta .* zeta - 2 * xi .* zeta .* zeta + 2 .* eta .* zeta .* zeta + 4 * xi .* eta .* zeta - zeta .* zeta - 2 * xi .* zeta + 2 * eta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 7) = (eta .* zeta .* (2 * xi + 1) .* (zeta + 1)) / 8 + (zeta .* (2 * xi + 1) .* (eta + 1) .* (zeta + 1)) / 8;
            d2N_xi_xi_k_I(1, 2, :, 8) = (1 / 8) * (4 * xi .* eta .* zeta .* zeta + 2 * xi .* zeta .* zeta - 2 .* eta .* zeta .* zeta + 4 * xi .* eta .* zeta - zeta .* zeta + 2 * xi .* zeta - 2 * eta .* zeta - zeta);
            d2N_xi_xi_k_I(1, 2, :, 9) = 0.25 * (-4 * xi .* eta .* zeta .* zeta + 2 .* xi .* zeta .* zeta + 2 * eta .* zeta .* zeta - zeta .* zeta + 4 * xi .* eta - 2 * xi - 2 * eta + 1);
            d2N_xi_xi_k_I(1, 2, :, 10) = -(eta .* (2 * xi + 1) .* (zeta.^2 - 1)) / 4 - ((2 * xi + 1) .* (zeta.^2 - 1) .* (eta - 1)) / 4;
            d2N_xi_xi_k_I(1, 2, :, 11) = 0.25 * (-4 * xi .* eta .* zeta .* zeta - 2 .* xi .* zeta .* zeta - 2 * eta .* zeta .* zeta - zeta .* zeta + 4 * xi .* eta + 2 * xi + 2 * eta + 1);
            d2N_xi_xi_k_I(1, 2, :, 12) = 0.25 * (-4 * xi .* eta .* zeta .* zeta - 2 .* xi .* zeta .* zeta + 2 * eta .* zeta .* zeta + zeta .* zeta + 4 * xi .* eta + 2 * xi - 2 * eta - 1);
            d2N_xi_xi_k_I(1, 2, :, 13) = 0.5 * (-2 * xi .* eta .* zeta .* zeta + xi .* zeta .* zeta + 2 * xi .* eta .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 14) = 0.5 * (-2 * xi .* eta .* zeta .* zeta - eta .* zeta .* zeta + 2 * xi .* eta .* zeta + eta .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 15) = 0.5 * (-2 * xi .* eta .* zeta .* zeta - xi .* zeta .* zeta + 2 * xi .* eta .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 16) = 0.5 * (-2 * xi .* eta .* zeta .* zeta + eta .* zeta .* zeta + 2 * xi .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 17) = 2 * xi .* eta .* zeta .* zeta - 2 .* xi .* eta .* zeta;
            d2N_xi_xi_k_I(1, 2, :, 18) = 0.5 * (-2 * xi .* eta .* zeta .* zeta + xi .* zeta .* zeta - 2 * xi .* eta .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 19) = 0.5 * (-2 * xi .* eta .* zeta .* zeta - eta .* zeta .* zeta - 2 * xi .* eta .* zeta - eta .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 20) = 0.5 * (-2 * xi .* eta .* zeta .* zeta - xi .* zeta .* zeta - 2 * xi .* eta .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(1, 2, :, 21) = -(eta .* zeta .* (2 * xi - 1) .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(1, 2, :, 22) = 2 * xi .* eta .* zeta .* zeta + 2 .* xi .* eta .* zeta;
            d2N_xi_xi_k_I(1, 2, :, 23) = 2 * xi .* eta .* zeta .* zeta - xi .* zeta .* zeta - 2 * xi .* eta + xi;
            d2N_xi_xi_k_I(1, 2, :, 24) = 2 * xi .* eta .* zeta .* zeta + eta .* zeta .* zeta - 2 * xi .* eta - eta;
            d2N_xi_xi_k_I(1, 2, :, 25) = 2 * xi .* eta .* zeta .* zeta + xi .* zeta .* zeta - 2 * xi .* eta - xi;
            d2N_xi_xi_k_I(1, 2, :, 26) = 2 * xi .* eta .* zeta .* zeta - eta .* zeta .* zeta - 2 * xi .* eta + eta;
            d2N_xi_xi_k_I(1, 2, :, 27) = -4 * xi .* eta .* zeta .* zeta + 4 * xi .* eta;
            d2N_xi_xi_k_I(1, 3, :, 1) = (1 / 8) * (4 * xi .* eta .* eta .* zeta - 2 * xi .* eta .* eta - 2 .* eta .* eta .* zeta - 4 * xi .* eta .* zeta + eta .* eta + 2 * xi .* eta + 2 * eta .* zeta - eta);
            d2N_xi_xi_k_I(1, 3, :, 2) = (1 / 8) * (4 * xi .* eta .* eta .* zeta - 2 * xi .* eta .* eta + 2 .* eta .* eta .* zeta - 4 * xi .* eta .* zeta - eta .* eta + 2 * xi .* eta - 2 * eta .* zeta + eta);
            d2N_xi_xi_k_I(1, 3, :, 3) = (eta .* zeta .* (2 * xi + 1) .* (eta + 1)) / 8 + (eta .* (2 * xi + 1) .* (eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 3, :, 4) = (eta .* zeta .* (2 * xi - 1) .* (eta + 1)) / 8 + (eta .* (2 * xi - 1) .* (eta + 1) .* (zeta - 1)) / 8;
            d2N_xi_xi_k_I(1, 3, :, 5) = (1 / 8) * (4 * xi .* eta .* eta .* zeta + 2 * xi .* eta .* eta - 2 .* eta .* eta .* zeta - 4 * xi .* eta .* zeta - eta .* eta - 2 * xi .* eta + 2 * eta .* zeta + eta);
            d2N_xi_xi_k_I(1, 3, :, 6) = (1 / 8) * (4 * xi .* eta .* eta .* zeta + 2 * xi .* eta .* eta + 2 .* eta .* eta .* zeta - 4 * xi .* eta .* zeta + eta .* eta - 2 * xi .* eta - 2 * eta .* zeta - eta);
            d2N_xi_xi_k_I(1, 3, :, 7) = (1 / 8) * (4 * xi .* eta .* eta .* zeta + 2 * xi .* eta .* eta + 2 .* eta .* eta .* zeta + 4 * xi .* eta .* zeta + eta .* eta + 2 * xi .* eta + 2 * eta .* zeta + eta);
            d2N_xi_xi_k_I(1, 3, :, 8) = (1 / 8) * (4 * xi .* eta .* eta .* zeta + 2 * xi .* eta .* eta - 2 .* eta .* eta .* zeta + 4 * xi .* eta .* zeta - eta .* eta + 2 * xi .* eta - 2 * eta .* zeta - eta);
            d2N_xi_xi_k_I(1, 3, :, 9) = -(eta .* zeta .* (2 * xi - 1) .* (eta - 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 10) = -(eta .* zeta .* (2 * xi + 1) .* (eta - 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 11) = -(eta .* zeta .* (2 * xi + 1) .* (eta + 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 12) = -(eta .* zeta .* (2 * xi - 1) .* (eta + 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 13) = 0.5 * (-2 * xi .* eta .* eta .* zeta + xi .* eta .* eta + 2 * xi .* eta .* zeta - xi .* eta);
            d2N_xi_xi_k_I(1, 3, :, 14) = 0.25 * (-4 * xi .* eta .* eta .* zeta + 2 .* xi .* eta .* eta - 2 * eta .* eta .* zeta + eta .* eta + 4 * xi .* zeta - 2 * xi + 2 * zeta - 1);
            d2N_xi_xi_k_I(1, 3, :, 15) = -(eta .* xi .* zeta .* (eta + 1)) / 2 - (eta .* xi .* (eta + 1) .* (zeta - 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 16) = -(zeta .* (eta.^2 - 1) .* (2 * xi - 1)) / 4 - ((eta.^2 - 1) .* (2 * xi - 1) .* (zeta - 1)) / 4;
            d2N_xi_xi_k_I(1, 3, :, 17) = 2 .* xi .* eta .* eta .* zeta - xi .* eta .* eta - 2 * xi .* zeta + xi;
            d2N_xi_xi_k_I(1, 3, :, 18) = -(eta .* xi .* zeta .* (eta - 1)) / 2 - (eta .* xi .* (eta - 1) .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(1, 3, :, 19) = 0.25 * (-4 * xi .* eta .* eta .* zeta - 2 .* xi .* eta .* eta - 2 * eta .* eta .* zeta - eta .* eta + 4 * xi .* zeta + 2 * xi + 2 * zeta + 1);
            d2N_xi_xi_k_I(1, 3, :, 20) = 0.5 * (-2 * xi .* eta .* eta .* zeta - xi .* eta .* eta - 2 * xi .* eta .* zeta - xi .* eta);
            d2N_xi_xi_k_I(1, 3, :, 21) = -(zeta .* (eta.^2 - 1) .* (2 * xi - 1)) / 4 - ((eta.^2 - 1) .* (2 * xi - 1) .* (zeta + 1)) / 4;
            d2N_xi_xi_k_I(1, 3, :, 22) = 2 .* xi .* eta .* eta .* zeta + xi .* eta .* eta - 2 * xi .* zeta - xi;
            d2N_xi_xi_k_I(1, 3, :, 23) = 2 .* xi .* eta .* eta .* zeta - 2 * xi .* eta .* zeta;
            d2N_xi_xi_k_I(1, 3, :, 24) = 2 .* xi .* eta .* eta .* zeta + eta .* eta .* zeta - 2 * xi .* zeta - zeta;
            d2N_xi_xi_k_I(1, 3, :, 25) = 2 .* xi .* eta .* eta .* zeta + 2 * xi .* eta .* zeta;
            d2N_xi_xi_k_I(1, 3, :, 26) = 2 .* xi .* eta .* eta .* zeta - eta .* eta .* zeta - 2 * xi .* zeta + zeta;
            d2N_xi_xi_k_I(1, 3, :, 27) = -4 * xi .* eta .* eta .* zeta + 4 * xi .* zeta;
            % second derivative wrt eta
            d2N_xi_xi_k_I(2, 1, :, 1) = d2N_xi_xi_k_I(1, 2, :, 1);
            d2N_xi_xi_k_I(2, 1, :, 2) = d2N_xi_xi_k_I(1, 2, :, 2);
            d2N_xi_xi_k_I(2, 1, :, 3) = d2N_xi_xi_k_I(1, 2, :, 3);
            d2N_xi_xi_k_I(2, 1, :, 4) = d2N_xi_xi_k_I(1, 2, :, 4);
            d2N_xi_xi_k_I(2, 1, :, 5) = d2N_xi_xi_k_I(1, 2, :, 5);
            d2N_xi_xi_k_I(2, 1, :, 6) = d2N_xi_xi_k_I(1, 2, :, 6);
            d2N_xi_xi_k_I(2, 1, :, 7) = d2N_xi_xi_k_I(1, 2, :, 7);
            d2N_xi_xi_k_I(2, 1, :, 8) = d2N_xi_xi_k_I(1, 2, :, 8);
            d2N_xi_xi_k_I(2, 1, :, 9) = d2N_xi_xi_k_I(1, 2, :, 9);
            d2N_xi_xi_k_I(2, 1, :, 10) = d2N_xi_xi_k_I(1, 2, :, 10);
            d2N_xi_xi_k_I(2, 1, :, 11) = d2N_xi_xi_k_I(1, 2, :, 11);
            d2N_xi_xi_k_I(2, 1, :, 12) = d2N_xi_xi_k_I(1, 2, :, 12);
            d2N_xi_xi_k_I(2, 1, :, 13) = d2N_xi_xi_k_I(1, 2, :, 13);
            d2N_xi_xi_k_I(2, 1, :, 14) = d2N_xi_xi_k_I(1, 2, :, 14);
            d2N_xi_xi_k_I(2, 1, :, 15) = d2N_xi_xi_k_I(1, 2, :, 15);
            d2N_xi_xi_k_I(2, 1, :, 16) = d2N_xi_xi_k_I(1, 2, :, 16);
            d2N_xi_xi_k_I(2, 1, :, 17) = d2N_xi_xi_k_I(1, 2, :, 17);
            d2N_xi_xi_k_I(2, 1, :, 18) = d2N_xi_xi_k_I(1, 2, :, 18);
            d2N_xi_xi_k_I(2, 1, :, 19) = d2N_xi_xi_k_I(1, 2, :, 19);
            d2N_xi_xi_k_I(2, 1, :, 20) = d2N_xi_xi_k_I(1, 2, :, 20);
            d2N_xi_xi_k_I(2, 1, :, 21) = d2N_xi_xi_k_I(1, 2, :, 21);
            d2N_xi_xi_k_I(2, 1, :, 22) = d2N_xi_xi_k_I(1, 2, :, 22);
            d2N_xi_xi_k_I(2, 1, :, 23) = d2N_xi_xi_k_I(1, 2, :, 23);
            d2N_xi_xi_k_I(2, 1, :, 24) = d2N_xi_xi_k_I(1, 2, :, 24);
            d2N_xi_xi_k_I(2, 1, :, 25) = d2N_xi_xi_k_I(1, 2, :, 25);
            d2N_xi_xi_k_I(2, 1, :, 26) = d2N_xi_xi_k_I(1, 2, :, 26);
            d2N_xi_xi_k_I(2, 1, :, 27) = d2N_xi_xi_k_I(1, 2, :, 27);
            d2N_xi_xi_k_I(2, 2, :, 1) = 0.25 * (xi .* xi .* zeta .* zeta - xi .* zeta .* zeta - xi .* xi .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 2) = 0.25 * (xi .* xi .* zeta .* zeta + xi .* zeta .* zeta - xi .* xi .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 3) = 0.25 * (xi .* xi .* zeta .* zeta + xi .* zeta .* zeta - xi .* xi .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 4) = 0.25 * (xi .* xi .* zeta .* zeta - xi .* zeta .* zeta - xi .* xi .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 5) = 0.25 * (xi .* xi .* zeta .* zeta - xi .* zeta .* zeta + xi .* xi .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 6) = 0.25 * (xi .* xi .* zeta .* zeta + xi .* zeta .* zeta + xi .* xi .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 7) = 0.25 * (xi .* xi .* zeta .* zeta + xi .* zeta .* zeta + xi .* xi .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 8) = 0.25 * (xi .* xi .* zeta .* zeta - xi .* zeta .* zeta + xi .* xi .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 9) = 0.5 * (-xi .* xi .* zeta .* zeta + xi .* zeta .* zeta + xi .* xi - xi);
            d2N_xi_xi_k_I(2, 2, :, 10) = -(xi .* (zeta.^2 - 1) .* (xi + 1)) / 2;
            d2N_xi_xi_k_I(2, 2, :, 11) = 0.5 * (-xi .* xi .* zeta .* zeta - xi .* zeta .* zeta + xi .* xi + xi);
            d2N_xi_xi_k_I(2, 2, :, 12) = 0.5 * (-xi .* xi .* zeta .* zeta + xi .* zeta .* zeta + xi .* xi - xi);
            d2N_xi_xi_k_I(2, 2, :, 13) = 0.5 * (-xi .* xi .* zeta .* zeta + xi .* xi .* zeta + zeta .* zeta - zeta);
            d2N_xi_xi_k_I(2, 2, :, 14) = -(xi .* zeta .* (xi + 1) .* (zeta - 1)) / 2;
            d2N_xi_xi_k_I(2, 2, :, 15) = 0.5 * (-xi .* xi .* zeta .* zeta + xi .* xi .* zeta + zeta .* zeta - zeta);
            d2N_xi_xi_k_I(2, 2, :, 16) = -(xi .* zeta .* (xi - 1) .* (zeta - 1)) / 2;
            d2N_xi_xi_k_I(2, 2, :, 17) = xi .* xi .* zeta .* zeta - xi .* xi .* zeta - zeta .* zeta + zeta;
            d2N_xi_xi_k_I(2, 2, :, 18) = 0.5 * (-xi .* xi .* zeta .* zeta - xi .* xi .* zeta + zeta .* zeta + zeta);
            d2N_xi_xi_k_I(2, 2, :, 19) = -(xi .* zeta .* (xi + 1) .* (zeta + 1)) / 2;
            d2N_xi_xi_k_I(2, 2, :, 20) = 0.5 * (-xi .* xi .* zeta .* zeta - xi .* xi .* zeta + zeta .* zeta + zeta);
            d2N_xi_xi_k_I(2, 2, :, 21) = 0.5 * (-xi .* xi .* zeta .* zeta + xi .* zeta .* zeta - xi .* xi .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 2, :, 22) = xi .* xi .* zeta .* zeta + xi .* xi .* zeta - zeta .* zeta - zeta;
            d2N_xi_xi_k_I(2, 2, :, 23) = xi .* xi .* zeta .* zeta - xi .* xi - zeta .* zeta + 1;
            d2N_xi_xi_k_I(2, 2, :, 24) = xi .* (zeta.^2 - 1) .* (xi + 1);
            d2N_xi_xi_k_I(2, 2, :, 25) = xi .* xi .* zeta .* zeta - xi .* xi - zeta .* zeta + 1;
            d2N_xi_xi_k_I(2, 2, :, 26) = xi .* xi .* zeta .* zeta - xi .* zeta .* zeta - xi .* xi + xi;
            d2N_xi_xi_k_I(2, 2, :, 27) = -2 * xi .* xi .* zeta .* zeta + 2 * xi .* xi + 2 * zeta .* zeta - 2;
            d2N_xi_xi_k_I(2, 3, :, 1) = (1 / 8) * (4 * xi .* xi .* eta .* zeta - 2 * xi .* xi .* eta - 2 * xi .* xi .* zeta - 4 * xi .* eta .* zeta + xi .* xi + 2 * xi .* eta + 2 * xi .* zeta - xi);
            d2N_xi_xi_k_I(2, 3, :, 2) = (1 / 8) * (4 * xi .* xi .* eta .* zeta - 2 * xi .* xi .* eta - 2 * xi .* xi .* zeta + 4 * xi .* eta .* zeta + xi .* xi - 2 * xi .* eta - 2 * xi .* zeta + xi);
            d2N_xi_xi_k_I(2, 3, :, 3) = (1 / 8) * (4 * xi .* xi .* eta .* zeta - 2 * xi .* xi .* eta + 2 * xi .* xi .* zeta + 4 * xi .* eta .* zeta - xi .* xi - 2 * xi .* eta + 2 * xi .* zeta - xi);
            d2N_xi_xi_k_I(2, 3, :, 4) = (1 / 8) * (4 * xi .* xi .* eta .* zeta - 2 * xi .* xi .* eta + 2 * xi .* xi .* zeta - 4 * xi .* eta .* zeta - xi .* xi + 2 * xi .* eta - 2 * xi .* zeta + xi);
            d2N_xi_xi_k_I(2, 3, :, 5) = (1 / 8) * (4 * xi .* xi .* eta .* zeta + 2 * xi .* xi .* eta - 2 * xi .* xi .* zeta - 4 * xi .* eta .* zeta - xi .* xi - 2 * xi .* eta + 2 * xi .* zeta + xi);
            d2N_xi_xi_k_I(2, 3, :, 6) = (xi .* zeta .* (2 * eta - 1) .* (xi + 1)) / 8 + (xi .* (2 * eta - 1) .* (xi + 1) .* (zeta + 1)) / 8;
            d2N_xi_xi_k_I(2, 3, :, 7) = (1 / 8) * (4 * xi .* xi .* eta .* zeta + 2 * xi .* xi .* eta + 2 * xi .* xi .* zeta + 4 * xi .* eta .* zeta + xi .* xi + 2 * xi .* eta + 2 * xi .* zeta + xi);
            d2N_xi_xi_k_I(2, 3, :, 8) = (1 / 8) * (4 * xi .* xi .* eta .* zeta + 2 * xi .* xi .* eta + 2 * xi .* xi .* zeta - 4 * xi .* eta .* zeta + xi .* xi - 2 * xi .* eta - 2 * xi .* zeta - xi);
            d2N_xi_xi_k_I(2, 3, :, 9) = 0.5 * (-2 * xi .* xi .* eta .* zeta + xi .* xi .* zeta + 2 * xi .* eta .* zeta - xi .* zeta);
            d2N_xi_xi_k_I(2, 3, :, 10) = 0.5 * (-2 * xi .* xi .* eta .* zeta + xi .* xi .* zeta - 2 * xi .* eta .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 3, :, 11) = -(xi .* zeta .* (2 * eta + 1) .* (xi + 1)) / 2;
            d2N_xi_xi_k_I(2, 3, :, 12) = 0.5 * (-2 * xi .* xi .* eta .* zeta - xi .* xi .* zeta + 2 * xi .* eta .* zeta + xi .* zeta);
            d2N_xi_xi_k_I(2, 3, :, 13) = -(zeta .* (2 * eta - 1) .* (xi.^2 - 1)) / 4 - ((2 * eta - 1) .* (xi.^2 - 1) .* (zeta - 1)) / 4;
            d2N_xi_xi_k_I(2, 3, :, 14) = 0.5 * (-2 * xi .* xi .* eta .* zeta + xi .* xi .* eta - 2 * xi .* eta .* zeta + xi .* eta);
            d2N_xi_xi_k_I(2, 3, :, 15) = -(zeta .* (2 * eta + 1) .* (xi.^2 - 1)) / 4 - ((2 * eta + 1) .* (xi.^2 - 1) .* (zeta - 1)) / 4;
            d2N_xi_xi_k_I(2, 3, :, 16) = 0.5 * (-2 * xi .* xi .* eta .* zeta + xi .* xi .* eta + 2 * xi .* eta .* zeta - xi .* eta);
            d2N_xi_xi_k_I(2, 3, :, 17) = 2 * xi .* xi .* eta .* zeta - xi .* xi .* eta - 2 * eta .* zeta + eta;
            d2N_xi_xi_k_I(2, 3, :, 18) = -(zeta .* (2 * eta - 1) .* (xi.^2 - 1)) / 4 - ((2 * eta - 1) .* (xi.^2 - 1) .* (zeta + 1)) / 4;
            d2N_xi_xi_k_I(2, 3, :, 19) = 0.5 * (-2 * xi .* xi .* eta .* zeta - xi .* xi .* eta - 2 * xi .* eta .* zeta - xi .* eta);
            d2N_xi_xi_k_I(2, 3, :, 20) = -(zeta .* (2 * eta + 1) .* (xi.^2 - 1)) / 4 - ((2 * eta + 1) .* (xi.^2 - 1) .* (zeta + 1)) / 4;
            d2N_xi_xi_k_I(2, 3, :, 21) = 0.5 * (-2 * xi .* xi .* eta .* zeta - xi .* xi .* eta + 2 * xi .* eta .* zeta + xi .* eta);
            d2N_xi_xi_k_I(2, 3, :, 22) = 2 * xi .* xi .* eta .* zeta + xi .* xi .* eta - 2 * eta .* zeta - eta;
            d2N_xi_xi_k_I(2, 3, :, 23) = 2 * xi .* xi .* eta .* zeta - xi .* xi .* zeta - 2 * eta .* zeta + zeta;
            d2N_xi_xi_k_I(2, 3, :, 24) = 2 * xi .* xi .* eta .* zeta + 2 * xi .* eta .* zeta;
            d2N_xi_xi_k_I(2, 3, :, 25) = 2 * xi .* xi .* eta .* zeta + xi .* xi .* zeta - 2 * eta .* zeta - zeta;
            d2N_xi_xi_k_I(2, 3, :, 26) = 2 * xi .* xi .* eta .* zeta - 2 * xi .* eta .* zeta;
            d2N_xi_xi_k_I(2, 3, :, 27) = -4 .* xi .* xi .* eta .* zeta + 4 * eta .* zeta;

            % second derivative wrt zeta
            d2N_xi_xi_k_I(3, 1, :, 1) = d2N_xi_xi_k_I(1, 3, :, 1);
            d2N_xi_xi_k_I(3, 1, :, 2) = d2N_xi_xi_k_I(1, 3, :, 2);
            d2N_xi_xi_k_I(3, 1, :, 3) = d2N_xi_xi_k_I(1, 3, :, 3);
            d2N_xi_xi_k_I(3, 1, :, 4) = d2N_xi_xi_k_I(1, 3, :, 4);
            d2N_xi_xi_k_I(3, 1, :, 5) = d2N_xi_xi_k_I(1, 3, :, 5);
            d2N_xi_xi_k_I(3, 1, :, 6) = d2N_xi_xi_k_I(1, 3, :, 6);
            d2N_xi_xi_k_I(3, 1, :, 7) = d2N_xi_xi_k_I(1, 3, :, 7);
            d2N_xi_xi_k_I(3, 1, :, 8) = d2N_xi_xi_k_I(1, 3, :, 8);
            d2N_xi_xi_k_I(3, 1, :, 9) = d2N_xi_xi_k_I(1, 3, :, 9);
            d2N_xi_xi_k_I(3, 1, :, 10) = d2N_xi_xi_k_I(1, 3, :, 10);
            d2N_xi_xi_k_I(3, 1, :, 11) = d2N_xi_xi_k_I(1, 3, :, 11);
            d2N_xi_xi_k_I(3, 1, :, 12) = d2N_xi_xi_k_I(1, 3, :, 12);
            d2N_xi_xi_k_I(3, 1, :, 13) = d2N_xi_xi_k_I(1, 3, :, 13);
            d2N_xi_xi_k_I(3, 1, :, 14) = d2N_xi_xi_k_I(1, 3, :, 14);
            d2N_xi_xi_k_I(3, 1, :, 15) = d2N_xi_xi_k_I(1, 3, :, 15);
            d2N_xi_xi_k_I(3, 1, :, 16) = d2N_xi_xi_k_I(1, 3, :, 16);
            d2N_xi_xi_k_I(3, 1, :, 17) = d2N_xi_xi_k_I(1, 3, :, 17);
            d2N_xi_xi_k_I(3, 1, :, 18) = d2N_xi_xi_k_I(1, 3, :, 18);
            d2N_xi_xi_k_I(3, 1, :, 19) = d2N_xi_xi_k_I(1, 3, :, 19);
            d2N_xi_xi_k_I(3, 1, :, 20) = d2N_xi_xi_k_I(1, 3, :, 20);
            d2N_xi_xi_k_I(3, 1, :, 21) = d2N_xi_xi_k_I(1, 3, :, 21);
            d2N_xi_xi_k_I(3, 1, :, 22) = d2N_xi_xi_k_I(1, 3, :, 22);
            d2N_xi_xi_k_I(3, 1, :, 23) = d2N_xi_xi_k_I(1, 3, :, 23);
            d2N_xi_xi_k_I(3, 1, :, 24) = d2N_xi_xi_k_I(1, 3, :, 24);
            d2N_xi_xi_k_I(3, 1, :, 25) = d2N_xi_xi_k_I(1, 3, :, 25);
            d2N_xi_xi_k_I(3, 1, :, 26) = d2N_xi_xi_k_I(1, 3, :, 26);
            d2N_xi_xi_k_I(3, 1, :, 27) = d2N_xi_xi_k_I(1, 3, :, 27);
            d2N_xi_xi_k_I(3, 2, :, 1) = d2N_xi_xi_k_I(2, 3, :, 1);
            d2N_xi_xi_k_I(3, 2, :, 2) = d2N_xi_xi_k_I(2, 3, :, 2);
            d2N_xi_xi_k_I(3, 2, :, 3) = d2N_xi_xi_k_I(2, 3, :, 3);
            d2N_xi_xi_k_I(3, 2, :, 4) = d2N_xi_xi_k_I(2, 3, :, 4);
            d2N_xi_xi_k_I(3, 2, :, 5) = d2N_xi_xi_k_I(2, 3, :, 5);
            d2N_xi_xi_k_I(3, 2, :, 6) = d2N_xi_xi_k_I(2, 3, :, 6);
            d2N_xi_xi_k_I(3, 2, :, 7) = d2N_xi_xi_k_I(2, 3, :, 7);
            d2N_xi_xi_k_I(3, 2, :, 8) = d2N_xi_xi_k_I(2, 3, :, 8);
            d2N_xi_xi_k_I(3, 2, :, 9) = d2N_xi_xi_k_I(2, 3, :, 9);
            d2N_xi_xi_k_I(3, 2, :, 10) = d2N_xi_xi_k_I(2, 3, :, 10);
            d2N_xi_xi_k_I(3, 2, :, 11) = d2N_xi_xi_k_I(2, 3, :, 11);
            d2N_xi_xi_k_I(3, 2, :, 12) = d2N_xi_xi_k_I(2, 3, :, 12);
            d2N_xi_xi_k_I(3, 2, :, 13) = d2N_xi_xi_k_I(2, 3, :, 13);
            d2N_xi_xi_k_I(3, 2, :, 14) = d2N_xi_xi_k_I(2, 3, :, 14);
            d2N_xi_xi_k_I(3, 2, :, 15) = d2N_xi_xi_k_I(2, 3, :, 15);
            d2N_xi_xi_k_I(3, 2, :, 16) = d2N_xi_xi_k_I(2, 3, :, 16);
            d2N_xi_xi_k_I(3, 2, :, 17) = d2N_xi_xi_k_I(2, 3, :, 17);
            d2N_xi_xi_k_I(3, 2, :, 18) = d2N_xi_xi_k_I(2, 3, :, 18);
            d2N_xi_xi_k_I(3, 2, :, 19) = d2N_xi_xi_k_I(2, 3, :, 19);
            d2N_xi_xi_k_I(3, 2, :, 20) = d2N_xi_xi_k_I(2, 3, :, 20);
            d2N_xi_xi_k_I(3, 2, :, 21) = d2N_xi_xi_k_I(2, 3, :, 21);
            d2N_xi_xi_k_I(3, 2, :, 22) = d2N_xi_xi_k_I(2, 3, :, 22);
            d2N_xi_xi_k_I(3, 2, :, 23) = d2N_xi_xi_k_I(2, 3, :, 23);
            d2N_xi_xi_k_I(3, 2, :, 24) = d2N_xi_xi_k_I(2, 3, :, 24);
            d2N_xi_xi_k_I(3, 2, :, 25) = d2N_xi_xi_k_I(2, 3, :, 25);
            d2N_xi_xi_k_I(3, 2, :, 26) = d2N_xi_xi_k_I(2, 3, :, 26);
            d2N_xi_xi_k_I(3, 2, :, 27) = d2N_xi_xi_k_I(2, 3, :, 27);
            d2N_xi_xi_k_I(3, 3, :, 1) = 0.25 * (xi .* xi .* eta .* eta - xi .* eta .* eta - xi .* xi .* eta + xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 2) = 0.25 * (xi .* xi .* eta .* eta + xi .* eta .* eta - xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 3) = (eta .* xi .* (eta + 1) .* (xi + 1)) / 4;
            d2N_xi_xi_k_I(3, 3, :, 4) = 0.25 * (xi .* xi .* eta .* eta - xi .* eta .* eta + xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 5) = 0.25 * (xi .* xi .* eta .* eta - xi .* eta .* eta - xi .* xi .* eta + xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 6) = 0.25 * (xi .* xi .* eta .* eta + xi .* eta .* eta - xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 7) = 0.25 * (xi .* xi .* eta .* eta + xi .* eta .* eta + xi .* xi .* eta + xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 8) = 0.25 * (xi .* xi .* eta .* eta - xi .* eta .* eta + xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 9) = 0.5 * (-xi .* xi .* eta .* eta + xi .* eta .* eta + xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 10) = 0.5 * (-xi .* xi .* eta .* eta - xi .* eta .* eta + xi .* xi .* eta + xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 11) = 0.5 * (-xi .* xi .* eta .* eta - xi .* eta .* eta - xi .* xi .* eta - xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 12) = 0.5 * (-xi .* xi .* eta .* eta + xi .* eta .* eta - xi .* xi .* eta + xi .* eta);
            d2N_xi_xi_k_I(3, 3, :, 13) = 0.5 * (-xi .* xi .* eta .* eta + xi .* xi .* eta + eta .* eta - eta);
            d2N_xi_xi_k_I(3, 3, :, 14) = 0.5 * (-xi .* xi .* eta .* eta - xi .* eta .* eta + xi .* xi + xi);
            d2N_xi_xi_k_I(3, 3, :, 15) = 0.5 * (-xi .* xi .* eta .* eta - xi .* xi .* eta + eta .* eta + eta);
            d2N_xi_xi_k_I(3, 3, :, 16) = 0.5 * (-xi .* xi .* eta .* eta + xi .* eta .* eta + xi .* xi - xi);
            d2N_xi_xi_k_I(3, 3, :, 17) = xi .* xi .* eta .* eta - xi .* xi - eta .* eta + 1;
            d2N_xi_xi_k_I(3, 3, :, 18) = 0.5 * (-xi .* xi .* eta .* eta + xi .* xi .* eta + eta .* eta - eta);
            d2N_xi_xi_k_I(3, 3, :, 19) = 0.5 * (-xi .* xi .* eta .* eta - xi .* eta .* eta + xi .* xi + xi);
            d2N_xi_xi_k_I(3, 3, :, 20) = 0.5 * (-xi .* xi .* eta .* eta - xi .* xi .* eta + eta .* eta + eta);
            d2N_xi_xi_k_I(3, 3, :, 21) = 0.5 * (-xi .* xi .* eta .* eta + xi .* eta .* eta + xi .* xi - xi);
            d2N_xi_xi_k_I(3, 3, :, 22) = xi .* xi .* eta .* eta - xi .* xi - eta .* eta + 1;
            d2N_xi_xi_k_I(3, 3, :, 23) = xi .* xi .* eta .* eta - xi .* xi .* eta - eta .* eta + eta;
            d2N_xi_xi_k_I(3, 3, :, 24) = xi .* xi .* eta .* eta + xi .* eta .* eta - xi .* xi - xi;
            d2N_xi_xi_k_I(3, 3, :, 25) = xi .* xi .* eta .* eta + xi .* xi .* eta - eta .* eta - eta;
            d2N_xi_xi_k_I(3, 3, :, 26) = xi .* xi .* eta .* eta - xi .* eta .* eta - xi .* xi + xi;
            d2N_xi_xi_k_I(3, 3, :, 27) = -2 * xi .* xi .* eta .* eta + 2 * xi .* xi + 2 * eta .* eta - 2;
        else
            error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
        end
    end
else
    error('Not implemented yet for given dimension!');
end
end

function [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = automaticComputationLagrangeShapeFunctions(dimension, numberOfNodes, evaluationPoints)
order = nthroot(numberOfNodes, dimension) - 1;
numberOfNodes1D = order + 1;

% compute lagrange basis polynomials and their derivatives
l = cell(dimension, numberOfNodes1D);
dl = cell(dimension, numberOfNodes1D);
d2l = cell(dimension, numberOfNodes1D);
for i = 1:numberOfNodes1D
    [l{1, i}, dl{1, i}, d2l{1, i}] = lagrangeBasisPolynomial(numberOfNodes1D, i, evaluationPoints(1, :).');
    if dimension >= 2
        [l{2, i}, dl{2, i}, d2l{2, i}] = lagrangeBasisPolynomial(numberOfNodes1D, i, evaluationPoints(2, :).');
    end
    if dimension >= 3
        [l{3, i}, dl{3, i}, d2l{3, i}] = lagrangeBasisPolynomial(numberOfNodes1D, i, evaluationPoints(3, :).');
    end
end

% compute lagrange shape functions and their derivatives
indicesXi = repmat(1:numberOfNodes1D, 1, numberOfNodes1D^(dimension - 1));
if dimension == 1
    N_k_I = [l{1, indicesXi}];
    dN_xi_k_I(1, :, :) = [dl{1, indicesXi}];
    d2N_xi_xi_k_I(1, 1, :, :) = [d2l{1, indicesXi}];
elseif dimension == 2
    indicesEta = kron(1:numberOfNodes1D, ones(1, numberOfNodes1D^(dimension - 1)));

    N_k_I = [l{1, indicesXi}] .* [l{2, indicesEta}];
    dN_xi_k_I(1, :, :) = [dl{1, indicesXi}] .* [l{2, indicesEta}];
    dN_xi_k_I(2, :, :) = [l{1, indicesXi}] .* [dl{2, indicesEta}];
    d2N_xi_xi_k_I(1, 1, :, :) = [d2l{1, indicesXi}] .* [l{2, indicesEta}];
    d2N_xi_xi_k_I(1, 2, :, :) = [dl{1, indicesXi}] .* [dl{2, indicesEta}];
    d2N_xi_xi_k_I(2, 1, :, :) = d2N_xi_xi_k_I(1, 2, :, :);
    d2N_xi_xi_k_I(2, 2, :, :) = [l{1, indicesXi}] .* [d2l{2, indicesEta}];
elseif dimension == 3
    indicesEta = kron(repmat(1:numberOfNodes1D, 1, numberOfNodes1D^(dimension - 2)), ones(1, numberOfNodes1D^(dimension - 2)));
    indicesZeta = kron(1:numberOfNodes1D, ones(1, numberOfNodes1D^(dimension - 1)));

    N_k_I = [l{1, indicesXi}] .* [l{2, indicesEta}] .* [l{3, indicesZeta}];
    dN_xi_k_I(1, :, :) = [dl{1, indicesXi}] .* [l{2, indicesEta}] .* [l{3, indicesZeta}];
    dN_xi_k_I(2, :, :) = [l{1, indicesXi}] .* [dl{2, indicesEta}] .* [l{3, indicesZeta}];
    dN_xi_k_I(3, :, :) = [l{1, indicesXi}] .* [l{2, indicesEta}] .* [dl{3, indicesZeta}];
    d2N_xi_xi_k_I(1, 1, :, :) = [d2l{1, indicesXi}] .* [l{2, indicesEta}] .* [l{3, indicesZeta}];
    d2N_xi_xi_k_I(1, 2, :, :) = [dl{1, indicesXi}] .* [dl{2, indicesEta}] .* [l{3, indicesZeta}];
    d2N_xi_xi_k_I(1, 3, :, :) = [dl{1, indicesXi}] .* [l{2, indicesEta}] .* [dl{3, indicesZeta}];
    d2N_xi_xi_k_I(2, 1, :, :) = d2N_xi_xi_k_I(1, 2, :, :);
    d2N_xi_xi_k_I(2, 2, :, :) = [l{1, indicesXi}] .* [d2l{2, indicesEta}] .* [l{3, indicesZeta}];
    d2N_xi_xi_k_I(2, 3, :, :) = [l{1, indicesXi}] .* [dl{2, indicesEta}] .* [dl{3, indicesZeta}];
    d2N_xi_xi_k_I(3, 1, :, :) = d2N_xi_xi_k_I(1, 3, :, :);
    d2N_xi_xi_k_I(3, 2, :, :) = d2N_xi_xi_k_I(2, 3, :, :);
    d2N_xi_xi_k_I(3, 3, :, :) = [l{1, indicesXi}] .* [l{2, indicesEta}] .* [d2l{3, indicesZeta}];
end

% rearrange matrices according to edof
if dimension == 1
    [~, edof] = meshOneDimensional(2, 1, order);
elseif dimension == 2
    edof = meshEdofRectangle(1, 1, 1, order);
elseif dimension == 3
    [~, edof] = meshGeneratorCube(2, 2, 2, 1, 1, 1, order, false);
end
N_k_I = N_k_I(:, int32(edof));
dN_xi_k_I = dN_xi_k_I(:, :, int32(edof));
d2N_xi_xi_k_I = d2N_xi_xi_k_I(:, :, :, int32(edof));
end

function [l, dl, d2l] = lagrangeBasisPolynomial(numberOfNodes1D, j, evaluationPoint)
nodeCoordinates = linspace(-1, 1, numberOfNodes1D);

indices = true(size(nodeCoordinates));
indices(j) = false;

l = prod((evaluationPoint - nodeCoordinates(indices))./(nodeCoordinates(j) - nodeCoordinates(indices)), 2);
factor1 = sum(1./(evaluationPoint - nodeCoordinates(indices)), 2);
factor2 = sum(1./(evaluationPoint - nodeCoordinates(indices)).^2, 2);
dl = l .* factor1;
d2l = l .* (factor1.^2 - factor2);
end