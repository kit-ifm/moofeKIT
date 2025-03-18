function [rData, kData, elementEnergy, array] = displacementBatheDvorkin1985HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% DISPLACEMENTBATHEDVORKIN1985HOOKEENDPOINT Element routine of class plateClass.
%
% FORMULATION
% This is the classical Bubnov-Galerkin Bathe-Dvorkin plate element.
% It falls in the class of irreducible Assumed Natural Strain (ANS)
% formulations and is also known as MITC4 (Mixed Interpolation of Tensorial
% Components) element.
% Only suitable for materially and geometrically linear simulations.
%
% CALL
% displacementBatheDvorkin1985HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [1, 1] for residual data.
% kData: cell-array of size [1, 1] for tangent data.
% dofs: degrees of freedom (dofs) optionally manipulated data (numerical
%       tangent)
% array: structure for storage fe-data, for more information see
%        storageFEObject.initializeArrayStress
% stressTensor: structure for storage stress tensors (postprocessing), for
%               more information see storageFEObject.initializeArrayStress
% flagNumericalTangent: flag that indicates whether the function call
%                       happens during the computation of the numerical
%                       tangent or not.
%
% REFERENCE
% https://doi.org/10.1002/nme.1620210213
%
% SEE ALSO
% -
%
% CREATOR(S)
% Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

numberOfNodes = size(dN_xi_k_I, 3);
gaussPoints = shapeFunctionObject.gaussPoint;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% compute material matrices
% bending component
Eb = zeros(3, 3);
Eb(1, 1) = 1;
Eb(1, 2) = nu;
Eb(2, 1) = nu;
Eb(2, 2) = 1;
Eb(3, 3) = (1 - nu) / 2;
Eb = E * h^3 / (12 * (1 - nu^2)) * Eb;
% shear component
Es = eye(2, 2);
Es = 5 / 6 * E / (2 * (1 + nu)) * h * Es;

% initialize energy
elementEnergy.strainEnergy = 0;

% compute Jacobian
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

% shapeFunctions at collocation points
collocationPoints = [0, -1, 0, 1; 1, 0, -1, 0];
[NGamma_k_I, dNGamma_xi_k_I, ~] = computeLagrangeShapeFunction(2, 4, 4, collocationPoints);
JAllCollocationPoints = computeJacobianForAllGausspoints(ed, dNGamma_xi_k_I);
JA = extractJacobianForGausspoint(JAllCollocationPoints, 1, setupObject, dimension);
JB = extractJacobianForGausspoint(JAllCollocationPoints, 2, setupObject, dimension);
JC = extractJacobianForGausspoint(JAllCollocationPoints, 3, setupObject, dimension);
JD = extractJacobianForGausspoint(JAllCollocationPoints, 4, setupObject, dimension);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % B-Matrix
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);

    % shear component
    BsANS = zeros(2, 3*numberOfNodes);
    BsANS(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dNGamma_xi_k_I(1, 1, :) + (1 - gaussPoints(2, k)) * dNGamma_xi_k_I(1, 3, :));
    BsANS(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dNGamma_xi_k_I(2, 4, :) + (1 - gaussPoints(1, k)) * dNGamma_xi_k_I(2, 2, :));
    BsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * NGamma_k_I(3, :));
    BsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * NGamma_k_I(2, :));
    BsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * NGamma_k_I(3, :));
    BsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * NGamma_k_I(2, :));

    Bs = J.' \ BsANS;

    if ~computePostData
        % Tangent
        % bending component
        Kb = Bb' * Eb * Bb * detJ * gaussWeight(k);
        % shear component
        Ks = Bs' * Es * Bs * detJ * gaussWeight(k);

        KXX = KXX + Kb + Ks;
    else
        % stress at gausspoint
        kappa = Bb * qN1(:);
        gamma = Bs * qN1(:);

        m = Eb * kappa;
        q = Es * gamma;
        stressTensor.Cauchy = [m(1), m(3), q(1); m(3), m(2), q(2); q(1), q(2), 0];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

%% RESIDUAL
RX = RX + KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end