function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% DISPLACEMENTPETROVGALERKINHOOKEENDPOINT Element routine of class plateClass.
%
% FORMULATION
% This an irreducible Petrov-Galerkin plate element employing the standard
% lagrangian shape functions for the virtual displacements and metric shape
% functions for the displacement field.
% This formulation is only suitable for research, as it suffers strongly
% from transverse shear locking.
% Only suitable for materially and geometrically linear simulations.
%
% CALL
% displacementPetrovGalerkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

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

% compute Jacobian matrices
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
Es = 5 / 6 * E / (2 * (1 + nu)) * h * Es;

% initialize energy
elementEnergy.strainEnergy = 0;

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

Kb = zeros(size(KXX));
Ks = zeros(size(KXX));

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, numberOfNodes);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, ed, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, ed, J0, shapeFunctionObject);
[M_k_I, dM_xi_k_I] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, J] = computeAllJacobian(ed, ed, ed, dN_xi_k_I, k, setupObject);

    % derivatives with respect to the physical coordinates
    index = dimension * (k - 1) + 1:dimension * k;
    dM_X_k_I = J0' \ dM_xi_k_I(index, :);

    % B-Matrix
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);

    Lb = zeros(3, 3*numberOfNodes);
    Lb(1, 2:3:end) = dM_X_k_I(1, :);
    Lb(2, 3:3:end) = dM_X_k_I(2, :);
    Lb(3, 2:3:end) = dM_X_k_I(2, :);
    Lb(3, 3:3:end) = dM_X_k_I(1, :);

    % shear component
    Bs = zeros(2, 3*numberOfNodes);
    Bs(1, 1:3:end) = dN_X_I(1, :);
    Bs(2, 1:3:end) = dN_X_I(2, :);
    Bs(1, 2:3:end) = N_k_I(k, :);
    Bs(2, 3:3:end) = N_k_I(k, :);

    Ls = zeros(2, 3*numberOfNodes);
    Ls(1, 1:3:end) = dM_X_k_I(1, :);
    Ls(2, 1:3:end) = dM_X_k_I(2, :);
    Ls(1, 2:3:end) = M_k_I(k, :);
    Ls(2, 3:3:end) = M_k_I(k, :);

    if ~computePostData
        % Tangent
        % bending component
        Kb = Kb + Bb' * Eb * Lb * detJ * gaussWeight(k);
        % shear component
        Ks = Ks + Bs' * Es * Ls * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        kappa = Lb * qN1(:);
        gamma = Ls * qN1(:);

        m = Eb * kappa;
        q = Es * gamma;
        stressTensor.Cauchy = [m(1), m(3), q(1); m(3), m(2), q(2); q(1), q(2), 0];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

KXX = KXX + Kb + Ks;

%% RESIDUAL
RX = RX + KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end
