function [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I, d3N_xi_xi_xi_k_I] = computeCubicHermiteShapeFunctions(dimension, numberOfGausspoints, gaussPoints)
%COMPUTECUBICHERMITESHAPEFUNCTIONS Computes the ansatz functions with cubic hermite polynomials in 1D.
%
%   CALL
%   [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = computeCubicHermiteShapeFunctions(dimension, numberOfGausspoints, gaussPoints)
%   dimension: dimension (has to be 1!)
%   numberOfGausspoints: number of Gauss points
%   gaussPoint: Gauss points
%   N_k_I: Ansatz function
%   dN_xi_k_I: First derivatives of the ansatz functions w.r.t local
%              coordinate
%   d2N_xi_xi_k_I: Second derivatives of the ansatz function w.r.t local
%                  coordinates
%
%   NOTE:
%   This function is used in the element routines, where also the
%   transformation to the physical domain takes place.
%
%   CREATOR(S)
%   Philipp Kinon

% Only implemented for 1D structures (e.g. Bernoulli beam)
assert(dimension == 1,"Cubic Hermite Ansatz functions only for 1D implemented.")

% Initialize arrays
N_k_I = zeros(numberOfGausspoints, 4);
dN_xi_k_I = zeros(numberOfGausspoints, 4);
d2N_xi_xi_k_I = zeros(numberOfGausspoints, 4);

% Local coordinate evaluated at gausspoints
xi = gaussPoints;

% Ansatz polynomials
N_k_I(:, 1) = 1/4 * (2 - 3*xi + xi.^3);          % has value 1 and slope 0 at xi=-1 and value 0 and slope 0 at xi=+1
N_k_I(:, 2) = 1/4 * (1 - xi - xi.^2 + xi.^3);    % has value 0 and slope 1 at xi=-1 and value 0 and slope 0 at xi=+1
N_k_I(:, 3) = 1/4 * (2 + 3*xi - xi.^3);          % has value 0 and slope 0 at xi=-1 and value 1 and slope 0 at xi=+1
N_k_I(:, 4) = 1/4 * (-1 - xi + xi.^2 + xi.^3);   % has value 0 and slope 0 at xi=-1 and value 0 and slope 1 at xi=+1

% First derivatives
dN_xi_k_I(:, 1) = 1/4 * (-3 + 3*xi.^2);
dN_xi_k_I(:, 2) = 1/4 * (-1 - 2*xi +3*xi.^2);
dN_xi_k_I(:, 3) = 1/4 * ( 3 - 3*xi.^2);
dN_xi_k_I(:, 4) = 1/4 * (-1 + 2*xi +3*xi.^2);

% Second derivatives
d2N_xi_xi_k_I(:,1) = 3/2*xi;
d2N_xi_xi_k_I(:,2) = 3/2*xi - 1/2;
d2N_xi_xi_k_I(:,3) = -3/2*xi;
d2N_xi_xi_k_I(:,4) = 3/2*xi + 1/2;

% Third derivatives
d3N_xi_xi_xi_k_I(:,1) = 3/2*ones(1,length(xi));
d3N_xi_xi_xi_k_I(:,2) = 3/2*ones(1,length(xi));
d3N_xi_xi_xi_k_I(:,3) = -3/2*ones(1,length(xi));
d3N_xi_xi_xi_k_I(:,4) = 3/2*ones(1,length(xi));

end