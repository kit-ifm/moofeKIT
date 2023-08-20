function [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = computeLagrangeShapeFunctionSymbolic(dimension, numberOfNodes, numberOfGausspoints, gaussPoints)
%COMPUTELAGRANGESHAPEFUNCTIONSYMBOLIC Symbolic computation of Lagrange shape functions.
%   This function computes Lagrange shape functions and their derivatives
%   using MATLABs Symbolic Toolbox. This function can be used to check the
%   correctness of other functions with the same objective.

syms xi eta zeta;

% initialize cells / arrays
NSymbolic = cell(numberOfNodes, 1);
N_k_I = zeros(numberOfGausspoints, numberOfNodes);
dN_xi_k_I = zeros(dimension, numberOfGausspoints, numberOfNodes);
d2N_xi_xi_k_I = zeros(dimension, dimension, numberOfGausspoints, numberOfNodes);

% compute shape functions
if dimension == 1
    if numberOfNodes == 1
        NSymbolic{1} = 1;
    elseif numberOfNodes == 2
        NSymbolic{1} = 1 / 2 * (1 - xi);
        NSymbolic{2} = 1 / 2 * (1 + xi);
    elseif numberOfNodes == 3
        NSymbolic{1} = 1 / 2 * xi .* (xi - 1);
        NSymbolic{2} = 1 / 2 * xi .* (xi + 1);
        NSymbolic{3} = 1 - xi.^2;
    else
        error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
    end
elseif dimension == 2
    if numberOfNodes == 1
        NSymbolic{1} = 1;
    elseif numberOfNodes == 3
        NSymbolic{1} = 1 - xi - eta;
        NSymbolic{2} = xi;
        NSymbolic{3} = eta;
    elseif numberOfNodes == 4
        NSymbolic{1} = (1 - xi) .* (1 - eta) / 4;
        NSymbolic{2} = (1 + xi) .* (1 - eta) / 4;
        NSymbolic{3} = (1 + xi) .* (1 + eta) / 4;
        NSymbolic{4} = (1 - xi) .* (1 + eta) / 4;
    elseif numberOfNodes == 6
        lambda = 1 - xi - eta;
        NSymbolic{1} = lambda .* (2 * lambda - 1);
        NSymbolic{2} = xi .* (2 * xi - 1);
        NSymbolic{3} = eta .* (2 * eta - 1);
        NSymbolic{4} = 4 * xi .* lambda;
        NSymbolic{5} = 4 * xi .* eta;
        NSymbolic{6} = 4 * eta .* lambda;
    elseif numberOfNodes == 8
        NSymbolic{1} = (1 - xi) .* (1 - eta) .* (-xi - eta - 1) / 4;
        NSymbolic{2} = (1 + xi) .* (1 - eta) .* (xi - eta - 1) / 4;
        NSymbolic{3} = (1 + xi) .* (1 + eta) .* (xi + eta - 1) / 4;
        NSymbolic{4} = (1 - xi) .* (1 + eta) .* (-xi + eta - 1) / 4;
        NSymbolic{5} = (1 - xi.^2) .* (1 - eta) / 2;
        NSymbolic{6} = (1 + xi) .* (1 - eta.^2) / 2;
        NSymbolic{7} = (1 - xi.^2) .* (1 + eta) / 2;
        NSymbolic{8} = (1 - xi) .* (1 - eta.^2) / 2;
    elseif numberOfNodes == 9
        NSymbolic{1} = 1 / 4 * (xi.^2 - xi) .* (eta.^2 - eta);
        NSymbolic{2} = 1 / 4 * (xi.^2 + xi) .* (eta.^2 - eta);
        NSymbolic{3} = 1 / 4 * (xi.^2 + xi) .* (eta.^2 + eta);
        NSymbolic{4} = 1 / 4 * (xi.^2 - xi) .* (eta.^2 + eta);
        NSymbolic{5} = 1 / 2 * (1 - xi.^2) .* eta .* (eta - 1);
        NSymbolic{6} = 1 / 2 * xi .* (xi + 1) .* (1 - eta.^2);
        NSymbolic{7} = 1 / 2 * (1 - xi.^2) .* eta .* (eta + 1);
        NSymbolic{8} = 1 / 2 * xi .* (xi - 1) .* (1 - eta.^2);
        NSymbolic{9} = (1 - xi.^2) .* (1 - eta.^2);
    else
        error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
    end
elseif dimension == 3
    if numberOfNodes == 1
        NSymbolic{1} = 1;
    elseif numberOfNodes == 4
        NSymbolic{1} = 1 - xi - eta - zeta;
        NSymbolic{2} = xi;
        NSymbolic{3} = eta;
        NSymbolic{4} = zeta;
    elseif numberOfNodes == 8
        NSymbolic{1} = (1 - xi) .* (1 - eta) .* (1 - zeta) / 8;
        NSymbolic{2} = (1 + xi) .* (1 - eta) .* (1 - zeta) / 8;
        NSymbolic{3} = (1 + xi) .* (1 + eta) .* (1 - zeta) / 8;
        NSymbolic{4} = (1 - xi) .* (1 + eta) .* (1 - zeta) / 8;
        NSymbolic{5} = (1 - xi) .* (1 - eta) .* (1 + zeta) / 8;
        NSymbolic{6} = (1 + xi) .* (1 - eta) .* (1 + zeta) / 8;
        NSymbolic{7} = (1 + xi) .* (1 + eta) .* (1 + zeta) / 8;
        NSymbolic{8} = (1 - xi) .* (1 + eta) .* (1 + zeta) / 8;
    elseif numberOfNodes == 10
        lambda = 1 - xi - eta - zeta;
        NSymbolic{1} = lambda .* (2 * lambda - 1);
        NSymbolic{2} = xi .* (2 * xi - 1);
        NSymbolic{3} = eta .* (2 * eta - 1);
        NSymbolic{4} = zeta .* (2 * zeta - 1);
        NSymbolic{5} = 4 * xi .* lambda;
        NSymbolic{6} = 4 * xi .* eta;
        NSymbolic{7} = 4 * eta .* lambda;
        NSymbolic{8} = 4 * zeta .* lambda;
        NSymbolic{9} = 4 * xi .* zeta;
        NSymbolic{10} = 4 * eta .* zeta;
    elseif numberOfNodes == 20
        NSymbolic{1} = (1 - xi) .* (1 - eta) .* (1 - zeta) .* (-xi - eta - zeta - 2) / 8;
        NSymbolic{2} = (1 + xi) .* (1 - eta) .* (1 - zeta) .* (xi - eta - zeta - 2) / 8;
        NSymbolic{3} = (1 + xi) .* (1 + eta) .* (1 - zeta) .* (xi + eta - zeta - 2) / 8;
        NSymbolic{4} = (1 - xi) .* (1 + eta) .* (1 - zeta) .* (-xi + eta - zeta - 2) / 8;
        NSymbolic{5} = (1 - xi) .* (1 - eta) .* (1 + zeta) .* (-xi - eta + zeta - 2) / 8;
        NSymbolic{6} = (1 + xi) .* (1 - eta) .* (1 + zeta) .* (xi - eta + zeta - 2) / 8;
        NSymbolic{7} = (1 + xi) .* (1 + eta) .* (1 + zeta) .* (xi + eta + zeta - 2) / 8;
        NSymbolic{8} = (1 - xi) .* (1 + eta) .* (1 + zeta) .* (-xi + eta + zeta - 2) / 8;
        NSymbolic{9} = (1 - xi.^2) .* (1 - eta) .* (1 - zeta) / 4;
        NSymbolic{10} = (1 + xi) .* (1 - eta.^2) .* (1 - zeta) / 4;
        NSymbolic{11} = (1 - xi.^2) .* (1 + eta) .* (1 - zeta) / 4;
        NSymbolic{12} = (1 - xi) .* (1 - eta.^2) .* (1 - zeta) / 4;
        NSymbolic{13} = (1 - xi.^2) .* (1 - eta) .* (1 + zeta) / 4;
        NSymbolic{14} = (1 + xi) .* (1 - eta.^2) .* (1 + zeta) / 4;
        NSymbolic{15} = (1 - xi.^2) .* (1 + eta) .* (1 + zeta) / 4;
        NSymbolic{16} = (1 - xi) .* (1 - eta.^2) .* (1 + zeta) / 4;
        NSymbolic{17} = (1 - xi) .* (1 - eta) .* (1 - zeta.^2) / 4;
        NSymbolic{18} = (1 + xi) .* (1 - eta) .* (1 - zeta.^2) / 4;
        NSymbolic{19} = (1 + xi) .* (1 + eta) .* (1 - zeta.^2) / 4;
        NSymbolic{20} = (1 - xi) .* (1 + eta) .* (1 - zeta.^2) / 4;
    elseif numberOfNodes == 27
        NSymbolic{1} = xi .* eta .* zeta .* (xi - 1) .* (eta - 1) .* (zeta - 1) / 8;
        NSymbolic{2} = xi .* eta .* zeta .* (xi + 1) .* (eta - 1) .* (zeta - 1) / 8;
        NSymbolic{3} = xi .* eta .* zeta .* (xi + 1) .* (eta + 1) .* (zeta - 1) / 8;
        NSymbolic{4} = xi .* eta .* zeta .* (xi - 1) .* (eta + 1) .* (zeta - 1) / 8;
        NSymbolic{5} = xi .* eta .* zeta .* (xi - 1) .* (eta - 1) .* (zeta + 1) / 8;
        NSymbolic{6} = xi .* eta .* zeta .* (xi + 1) .* (eta - 1) .* (zeta + 1) / 8;
        NSymbolic{7} = xi .* eta .* zeta .* (xi + 1) .* (eta + 1) .* (zeta + 1) / 8;
        NSymbolic{8} = xi .* eta .* zeta .* (xi - 1) .* (eta + 1) .* (zeta + 1) / 8;
        NSymbolic{9} = xi .* eta .* (xi - 1) .* (eta - 1) .* (1 - zeta.^2) / 4;
        NSymbolic{10} = xi .* eta .* (xi + 1) .* (eta - 1) .* (1 - zeta.^2) / 4;
        NSymbolic{11} = xi .* eta .* (xi + 1) .* (eta + 1) .* (1 - zeta.^2) / 4;
        NSymbolic{12} = xi .* eta .* (xi - 1) .* (eta + 1) .* (1 - zeta.^2) / 4;
        NSymbolic{13} = eta .* zeta .* (1 - xi.^2) .* (eta - 1) .* (zeta - 1) / 4;
        NSymbolic{14} = xi .* zeta .* (xi + 1) .* (1 - eta.^2) .* (zeta - 1) / 4;
        NSymbolic{15} = eta .* zeta .* (1 - xi.^2) .* (eta + 1) .* (zeta - 1) / 4;
        NSymbolic{16} = xi .* zeta .* (xi - 1) .* (1 - eta.^2) .* (zeta - 1) / 4;
        NSymbolic{17} = zeta .* (1 - xi.^2) .* (1 - eta.^2) .* (zeta - 1) / 2;
        NSymbolic{18} = eta .* zeta .* (1 - xi.^2) .* (eta - 1) .* (zeta + 1) / 4;
        NSymbolic{19} = xi .* zeta .* (xi + 1) .* (1 - eta.^2) .* (zeta + 1) / 4;
        NSymbolic{20} = eta .* zeta .* (1 - xi.^2) .* (eta + 1) .* (zeta + 1) / 4;
        NSymbolic{21} = xi .* zeta .* (xi - 1) .* (1 - eta.^2) .* (zeta + 1) / 4;
        NSymbolic{22} = zeta .* (1 - xi.^2) .* (1 - eta.^2) .* (zeta + 1) / 2;
        NSymbolic{23} = eta .* (1 - xi.^2) .* (eta - 1) .* (1 - zeta.^2) / 2;
        NSymbolic{24} = xi .* (xi + 1) .* (1 - eta.^2) .* (1 - zeta.^2) / 2;
        NSymbolic{25} = eta .* (1 - xi.^2) .* (eta + 1) .* (1 - zeta.^2) / 2;
        NSymbolic{26} = xi .* (xi - 1) .* (1 - eta.^2) .* (1 - zeta.^2) / 2;
        NSymbolic{27} = (1 - xi.^2) .* (1 - eta.^2) .* (1 - zeta.^2);
    else
        error(['Not implemented yet for given numberOfNodes (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ')!']);
    end
else
    error('Not implemented yet for given dimension!');
end

% initialize values for substitution
if dimension == 1
    subsVars = {xi};
    subsValues = {gaussPoints(1, :)};
elseif dimension == 2
    subsVars = {xi, eta};
    subsValues = {gaussPoints(1, :), gaussPoints(2, :)};
elseif dimension == 3
    subsVars = {xi, eta, zeta};
    subsValues = {gaussPoints(1, :), gaussPoints(2, :), gaussPoints(3, :)};
else
    error('Not implemented yet for given dimension!');
end

% differentiation of shape functions and substitution
for ii = 1:numberOfNodes
    dNSymbolic_xi = diff(NSymbolic{ii}, xi);
    d2NSymbolic_xi_xi = diff(dNSymbolic_xi, xi);
    if dimension >= 2
        dNSymbolic_eta = diff(NSymbolic{ii}, eta);
        d2NSymbolic_xi_eta = diff(dNSymbolic_xi, eta);
        d2NSymbolic_eta_eta = diff(dNSymbolic_eta, eta);
    end
    if dimension >= 3
        dNSymbolic_zeta = diff(NSymbolic{ii}, zeta);
        d2NSymbolic_xi_zeta = diff(dNSymbolic_xi, zeta);
        d2NSymbolic_eta_zeta = diff(dNSymbolic_eta, zeta);
        d2NSymbolic_zeta_zeta = diff(dNSymbolic_zeta, zeta);
    end

    N_k_I(:, ii) = subs(NSymbolic{ii}, subsVars, subsValues);
    dN_xi_k_I(1, :, ii) = subs(dNSymbolic_xi, subsVars, subsValues);
    d2N_xi_xi_k_I(1, 1, :, ii) = subs(d2NSymbolic_xi_xi, subsVars, subsValues);
    if dimension >= 2
        dN_xi_k_I(2, :, ii) = subs(dNSymbolic_eta, subsVars, subsValues);
        d2N_xi_xi_k_I(1, 2, :, ii) = subs(d2NSymbolic_xi_eta, subsVars, subsValues);
        d2N_xi_xi_k_I(2, 1, :, ii) = d2N_xi_xi_k_I(1, 2, :, ii);
        d2N_xi_xi_k_I(2, 2, :, ii) = subs(d2NSymbolic_eta_eta, subsVars, subsValues);
    end
    if dimension >= 3
        dN_xi_k_I(3, :, ii) = subs(dNSymbolic_zeta, subsVars, subsValues);
        d2N_xi_xi_k_I(1, 3, :, ii) = subs(d2NSymbolic_xi_zeta, subsVars, subsValues);
        d2N_xi_xi_k_I(2, 3, :, ii) = subs(d2NSymbolic_eta_zeta, subsVars, subsValues);
        d2N_xi_xi_k_I(3, 1, :, ii) = d2N_xi_xi_k_I(1, 3, :, ii);
        d2N_xi_xi_k_I(3, 2, :, ii) = d2N_xi_xi_k_I(2, 3, :, ii);
        d2N_xi_xi_k_I(3, 3, :, ii) = subs(d2NSymbolic_zeta_zeta, subsVars, subsValues);
    end
end
end
