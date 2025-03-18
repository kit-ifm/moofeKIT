function [rData, kData, elementEnergy, array] = selectiveReducedIntegrationHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% SELECTIVEREDUCEDINTEGRATIONHOOKEENDPOINT Element routine of class plateClass.
%
% FORMULATION
% This is the classical Bubnov-Galerkin plate element employing selective
% reduced integration.
% Only suitable for materially and geometrically linear simulations.
%
% CALL
% selectiveReducedIntegrationHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% -
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
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;
N0_k_I = shapeFunctionObject.N0_k_I;

numberOfNodes = size(dN_xi_k_I, 3);
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

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

% compute Jacobian
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);
J0 = ed * dN0_xi_I';

%loop over all gauss points
for k = 1:numberOfGausspoints
    %approximation of geometry and deriviates according to x
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, 2);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k); %nodal operatormatix
    %B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    % B-MatrixBb
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);
    % shear component
    Bs = zeros(2, 3*numberOfNodes);
    Bs(:, 1:3:end) = J0.' \ dN0_xi_I;
    Bs(1, 2:3:end) = N0_k_I;
    Bs(2, 3:3:end) = N0_k_I;

    if ~computePostData
        %tangent and residuum
        Kb = Bb' * Eb * Bb * detJ * gaussWeight(k);
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
