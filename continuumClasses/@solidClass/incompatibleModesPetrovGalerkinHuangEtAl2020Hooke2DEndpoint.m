function [rData, kData, elementEnergy, array] = incompatibleModesPetrovGalerkinHuangEtAl2020Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% INCOMPATIBLEMODESPETROVGALERKINHUANGETAL2020HOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'incompatibleModes'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is IUQ4, where U denotes the unsymmetry of
% the stiffness matrix.
%
% CALL
% incompatibleModesPetrovGalerkinHuangEtAl2020Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data of
%        every field, here: (X, A)
% kData: cell-array of size [2, 2] for
%        tangent data of every field, here: (X, A)
% dofs: degrees of freedom (dofs)
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
% https://doi.org/10.1002/nme.6363
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

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

edof = meshObject.edof;

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;
switch lower(strtok(materialObject.name, 'Hooke'))
    case 'esz'
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end

% aquire the nodal values of the variables for the current element
X = obj.qR(edof(e, :), 1:dimension).';
x = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = x(:) - X(:);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);
J0 = X * dN0_xi_I';
detJ0 = det(J0);

% compute aVals and bVals
aVals = 1 / 4 * [-1, 1, 1, -1; 1, -1, 1, -1; -1, -1, 1, 1; 1, 1, 1, 1] * X(1, :).';
bVals = 1 / 4 * [-1, 1, 1, -1; 1, -1, 1, -1; -1, -1, 1, 1; 1, 1, 1, 1] * X(2, :).';

% compute additional shape functions
nodalPointsInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 9);
nodalPointsInPhysicalCoordinates = [aVals, bVals].' * [nodalPointsInLocalCoordinates(1, :); nodalPointsInLocalCoordinates(1, :) .* nodalPointsInLocalCoordinates(2, :); nodalPointsInLocalCoordinates(2, :); ones(1, 9)];
nodalPointsInSkewCoordinates = J0 \ (nodalPointsInPhysicalCoordinates - [aVals(4); bVals(4)]);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPoints, X, J0, shapeFunctionObject);
[~, dMBarr] = computeMetricShapeFunctions(obj, dimension, nodalPointsInSkewCoordinates, gaussPointsInSkewCoordinates);

dNBubble_xi_k_I = computeIncompatibleShapeFunctionsTestFunction(gaussPoints, numberOfGausspoints, detJ0, aVals, bVals);

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    dNBubble_X_I = computedN_X_I(dNBubble_xi_k_I, J, k);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMBarx = J0' \ dMBarr(indx, :);

    dMx = zeros(2, 4);
    dMx(:, 1) = dMBarx(:, 1) + 1 / 2 * (dMBarx(:, 5) + dMBarx(:, 8)) + 1 / 4 * dMBarx(:, 9);
    dMx(:, 2) = dMBarx(:, 2) + 1 / 2 * (dMBarx(:, 5) + dMBarx(:, 6)) + 1 / 4 * dMBarx(:, 9);
    dMx(:, 3) = dMBarx(:, 3) + 1 / 2 * (dMBarx(:, 6) + dMBarx(:, 7)) + 1 / 4 * dMBarx(:, 9);
    dMx(:, 4) = dMBarx(:, 4) + 1 / 2 * (dMBarx(:, 7) + dMBarx(:, 8)) + 1 / 4 * dMBarx(:, 9);

    dMBubblex = dMBarx(:, 5:9);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dN_X_I);
    L = BMatrix(dMx);
    G = BMatrix(dNBubble_X_I);
    H = BMatrix(dMBubblex);

    % strain tensor
    epsilonVoigt = L * uN1 + H * alphaN1e;

    if ~computePostData
        % TANGENT
        KDD = KDD + (B' * C * L) * detJ * gaussWeight(k);
        KDA = KDA + (B' * C * H) * detJ * gaussWeight(k);
        KAD = KAD + (G' * C * L) * detJ * gaussWeight(k);
        KAA = KAA + (G' * C * H) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaVoigt = C * epsilonVoigt;
        sigma = voigtToMatrix(sigmaVoigt, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
RD = KDD * uN1 + KDA * alphaN1e;
RA = KAD * uN1 + KAA * alphaN1e;

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RD;
    rData{2} = RA;
    kData{1, 1} = KDD;
    kData{1, 2} = KDA;
    kData{2, 1} = KAD;
    kData{2, 2} = KAA;
end
end

function dNBubble_xi_k_I = computeIncompatibleShapeFunctionsTestFunction(gaussPointsInLocalCoordinates, numberOfGausspoints, detJ0, aVals, bVals)
J1 = aVals(1) * bVals(2) - aVals(2) * bVals(1);
J2 = aVals(2) * bVals(3) - aVals(3) * bVals(2);

xi = gaussPointsInLocalCoordinates(1, :);
eta = gaussPointsInLocalCoordinates(2, :);

dNBubble_xi_k_I = zeros(2, numberOfGausspoints, 5);
dNBubble_xi_k_I(1, :, 1) = -xi .* eta .* (eta - 1) + J1 / (3 * detJ0) * ones(1, numberOfGausspoints);
dNBubble_xi_k_I(1, :, 2) = 1 / 2 * (2 * xi + 1) .* (1 - eta.^2) - (detJ0 + J1) / (3 * detJ0);
dNBubble_xi_k_I(1, :, 3) = -xi .* eta .* (eta + 1) + J1 / (3 * detJ0) * ones(1, numberOfGausspoints);
dNBubble_xi_k_I(1, :, 4) = 1 / 2 * (2 * xi - 1) .* (1 - eta.^2) - (-detJ0 + J1) / (3 * detJ0);
dNBubble_xi_k_I(1, :, 5) = -2 * xi .* (1 - eta.^2);

dNBubble_xi_k_I(2, :, 1) = 1 / 2 * (1 - xi.^2) .* (2 * eta - 1) - (-detJ0 + J2) / (3 * detJ0);
dNBubble_xi_k_I(2, :, 2) = -xi .* (xi + 1) .* eta + J2 / (3 * detJ0) * ones(1, numberOfGausspoints);
dNBubble_xi_k_I(2, :, 3) = 1 / 2 * (1 - xi.^2) .* (2 * eta + 1) - (detJ0 + J2) / (3 * detJ0);
dNBubble_xi_k_I(2, :, 4) = -xi .* (xi - 1) .* eta + J2 / (3 * detJ0) * ones(1, numberOfGausspoints);
dNBubble_xi_k_I(2, :, 5) = -2 * (1 - xi.^2) .* eta;
end