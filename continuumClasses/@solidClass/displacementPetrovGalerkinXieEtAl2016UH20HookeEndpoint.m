function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinXieEtAl2016UH20HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTPETROVGALERKINXIEETAL2016UH20HOOKEENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is UH20, where U denotes the unsymmetry of
% the stiffness matrix.
%
% CALL
% displacementPetrovGalerkinXieEtAl2016UH20Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% https://doi.org/10.1007/s10999-014-9289-3 (Xie, Sze, Zhou: Modified and Trefftz unsymmetric finite element models)
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

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;

numberOfNodes = size(N_k_I, 2);

dimension = obj.dimension;

% aquire material data
lambda = materialObject.lambda;
mu = materialObject.mu;
C = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
    lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
    lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
    0, 0, 0, mu, 0, 0; ...
    0, 0, 0, 0, mu, 0; ...
    0, 0, 0, 0, 0, mu];

% aquire the nodal values of the variables for the current element
X = obj.qR(edof(e, :), 1:dimension).';
x = dofs.edN1;
uN1 = x(:) - X(:);

% compute Jacobian matrices
J0 = X * dN0_xi_I';
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, numberOfNodes);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodalPointsInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, X, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, X, J0, shapeFunctionObject);
[~, dMr] = computeMetricShapeFunctions(obj, dimension, nodalPointsInSkewCoordinates, gaussPointsInSkewCoordinates);

% initialize residual
Re = rData{1};

% initialize tangent
Ke = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMx = J0' \ dMr(indx, :);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dN_X_I);
    L = BMatrix(dMx);
    
    if ~computePostData
        % TANGENT
        Ke = Ke + (B' * C * L) * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (L * uN1).' * C * (L * uN1) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaVoigt = C * L * uN1;
        sigma = voigtToMatrix(sigmaVoigt, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
Re = Re + Ke * uN1;

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end
