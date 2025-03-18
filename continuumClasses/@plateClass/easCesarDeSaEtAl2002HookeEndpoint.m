function [rData, kData, elementEnergy, array] = easCesarDeSaEtAl2002HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% EASCESARDESAETAL2002HOOKEENDPOINT Element routine of class plateClass.
%
% FORMULATION
% This is the classical Bubnov-Galerkin EAS plate element by Cesar de Sa et al.
% It is a mixed Enhanced Assumed Strain (EAS) formulation, where the
% transverse shear strains are enhanced to alleviate the transverse shear
% locking problem.
% Only suitable for regular meshes.
% Only suitable for materially and geometrically linear simulations.
% Does not pass the patch test!
%
% CALL
% easCesarDeSaEtAl2002HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data.
% kData: cell-array of size [2, 2] for tangent data.
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
% https://doi.org/10.1002/nme.360
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
mixedFEObject = obj.mixedFEObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;
ErAll = mixedFEObject.shapeFunctionObject.M;

numberOfNodes = size(dN_xi_k_I, 3);
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;
kappa = 5 / 6;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;
alphaN1 = dofs.edAlphaN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);
J0 = ed * dN0_xi_I';
detJ0 = det(J0);

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
Es = kappa * E / (2 * (1 + nu)) * h * Es;

% initialize energy
elementEnergy.strainEnergy = 0;

% initialize residual & tangent
RX = rData{1, 1};
RA = rData{2, 1};
KXX = kData{1, 1};
KXA = kData{1, 2};
KAX = kData{2, 1};
KAA = kData{2, 2};

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
    Bs = zeros(2, 3*numberOfNodes);
    Bs(1, 1:3:end) = dN_X_I(1, :);
    Bs(2, 1:3:end) = dN_X_I(2, :);
    Bs(1, 2:3:end) = N_k_I(k, :);
    Bs(2, 3:3:end) = N_k_I(k, :);

    % shape functions for the enhanced part of the strain field
    indx2 = 2 * (k - 1) + 1:2 * k;
    Er = ErAll(indx2, :);

    % approximation matrix
    G = detJ0 / detJ * (J0' \ Er);

    if ~computePostData
        % Tangent
        % bending component
        Kb = Bb' * Eb * Bb * detJ * gaussWeight(k);
        % shear component
        Ks = Bs' * Es * Bs * detJ * gaussWeight(k);

        KXX = KXX + Kb + Ks;
        KXA = KXA + Bs' * Es * G * detJ * gaussWeight(k);
        KAX = KAX + G' * Es * Bs * detJ * gaussWeight(k);
        KAA = KAA + G' * Es * G * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        kappa = Bb * qN1(:);
        gamma = Bs * qN1(:) + G * alphaN1(:);

        m = Eb * kappa;
        q = Es * gamma;
        stressTensor.Cauchy = [m(1), m(3), q(1); m(3), m(2), q(2); q(1), q(2), 0];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

%% RESIDUAL
RX = RX + KXX * qN1(:) + KXA * alphaN1(:);
RA = RA + KAX * qN1(:) + KAA * alphaN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    rData{2} = RA;
    kData{1, 1} = KXX;
    kData{1, 2} = KXA;
    kData{2, 1} = KAX;
    kData{2, 2} = KAA;
end
end