function [rData, kData, elementEnergy, array] = easPetrovGalerkinHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASPETROVGALERKINHOOKE2DENDPOINT Element routine of class solidThermoClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is Q1U/E4, where U denotes the unsymmetry of
% the stiffness matrix.
% This formulation extends the purely mechanical formulation to the case of
% steady-state thermoelasticity.
%
% CALL
% easPetrovGalerkinHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [3,1] for residual data of
%        every field, here: (X, T, alpha)
% kData: cell-array of size [3, 3] for
%        tangent data of every field, here: (X, T, alpha)
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
% solidClass/easPetrovGalerkinHooke2DEndpoint
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
ErAll = mixedFEObject.shapeFunctionObject.M;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;
alphaT = materialObject.alphaT;
kThermal = materialObject.k;
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
TRef = obj.qR(edof(e, :), dimension+1);
TN1 = dofs.thetaN1.';
uN1 = x(:) - X(:);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);
J0 = X * dN0_xi_I';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, X, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, X, J0, shapeFunctionObject);
[M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[~, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr);

% initialize tangent
KDD = kData{1, 1};
KDT = kData{1, 2};
KDA = kData{1, 3};
KTT = kData{2, 2};
KAD = kData{3, 1};
KAT = kData{3, 2};
KAA = kData{3, 3};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % shape functions for the enhanced part of the strain field
    indx2 = 3 * (k - 1) + 1:3 * k;
    Er = ErAll(indx2, :);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMx = J0' \ dMr(indx, :);
    dMTildex = J0' \ dMTilder(indx, :);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dN_X_I);
    L = BMatrix(dMx);
    G = 1 / detJ * (F0' \ Er);
    H = BMatrix(dMTildex);

    % strain tensor
    epsilonVoigt = L * uN1 + H * alphaN1e;
    epsilon0Voigt = alphaT * [1; 1; 0] * M_k_I(k, :) * (TN1-TRef);

    if ~computePostData
        % TANGENT
        KDD = KDD + (B' * C * L) * detJ * gaussWeight(k);
        KDA = KDA + (B' * C * H) * detJ * gaussWeight(k);
        KDT = KDT - (B' * C * alphaT * [1; 1; 0] * M_k_I(k, :)) * detJ * gaussWeight(k);
        KAD = KAD + (G' * C * L) * detJ * gaussWeight(k);
        KAA = KAA + (G' * C * H) * detJ * gaussWeight(k);
        KAT = KAT - (G' * C * alphaT * [1; 1; 0] * M_k_I(k, :)) * detJ * gaussWeight(k);
        KTT = KTT + (dN_X_I.' * kThermal * eye(2) * dMx) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaVoigt = C * (epsilonVoigt - epsilon0Voigt);
        sigma = voigtToMatrix(sigmaVoigt, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
RD = KDD * uN1 + KDA * alphaN1e + KDT * (TN1-TRef);
RT = KTT * (TN1-TRef);
RA = KAD * uN1 + KAA * alphaN1e + KAT * (TN1-TRef);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RD;
    rData{2} = RT;
    rData{3} = RA;
    kData{1, 1} = KDD;
    kData{1, 2} = KDT;
    kData{1, 3} = KDA;
    kData{2, 2} = KTT;
    kData{3, 1} = KAD;
    kData{3, 2} = KAT;
    kData{3, 3} = KAA;
end
end

function [MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr)
% Computes the shape functions for the trial function for the enhanced part
% of the strain field
MTilde = (gaussPointsInSkewCoordinates.^2).' - M_k_I * (nodesInSkewCoordinates.^2).';

dmTilder = zeros(8, 2);
dmTilder(1:2:end, 1) = 2 * gaussPointsInSkewCoordinates(1, :);
dmTilder(2:2:end, 2) = 2 * gaussPointsInSkewCoordinates(2, :);
dMTilder = dmTilder - dMr * (nodesInSkewCoordinates.^2).';
end