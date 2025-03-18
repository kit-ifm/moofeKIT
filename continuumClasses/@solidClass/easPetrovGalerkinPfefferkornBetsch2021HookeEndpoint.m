function [rData, kData, elementEnergy, array] = easPetrovGalerkinPfefferkornBetsch2021HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASPETROVGALERKINPFEFFERKORNBETSCH2021HOOKEENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is Q1U/E4 (2D) or H1U/E12, where U denotes
% the unsymmetry of the stiffness matrix.
%
% CALL
% easPetrovGalerkinPfefferkornBetsch2021HookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data of
%        every field, here: (X, alpha)
% kData: cell-array of size [2, 2] for
%        tangent data of every field, here: (X, alpha)
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
% https://doi.org/10.1002/nme.6817
%
% SEE ALSO
% easPetrovGalerkinPfefferkornBetsch2021Hooke2DEndpoint
% easHookeEndpoint
%
% CREATOR(S)
% Felix Zaehringer, Robin Pfefferkorn

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
easAnsatzFunctionsAll = mixedFEObject.shapeFunctionObject.M;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;

dimension = obj.dimension;

% aquire material data
if dimension == 2
    E = materialObject.E;
    nu = materialObject.nu;
    C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; ...
        nu, 1 - nu, 0; ...
        0, 0, (1 - 2 * nu) / 2];
else
    lambda = materialObject.lambda;
    mu = materialObject.mu;
    C = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
        lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
        lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
        0, 0, 0, mu, 0, 0; ...
        0, 0, 0, 0, mu, 0; ...
        0, 0, 0, 0, 0, mu];
end

% aquire the nodal values of the variables for the current element
X = obj.qR(edof(e, :), 1:dimension).';
xN = obj.qN(edof(e, :), 1:dimension).';
x = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = x(:) - X(:);

% compute Jacobian matrices
J0 = X * dN0_xi_I';
detJ0 = det(J0);
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 2^dimension);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, X, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, X, J0, shapeFunctionObject);
pianSumiharaAnsatzFunctionsAll = computePianSumiharaAnsatzFunctions(obj, dimension, 18, numberOfGausspoints, gaussPointsInSkewCoordinates);
[M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[~, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr);
ansatzFunctionTestFunctionEnhancedStrainAll = computeShapeFunctionsTestFunctionEnhancedStrain2(dimension, pianSumiharaAnsatzFunctionsAll, easAnsatzFunctionsAll, X, dN_xi_k_I, gaussWeight, nodesInLocalCoordinates);

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

    % shape functions for the enhanced part of the strain field
    indx2 = (3 * dimension - 3) * (k - 1) + 1:(3 * dimension - 3) * k;
    ansatzFunctionTestFunctionEnhancedStrain = ansatzFunctionTestFunctionEnhancedStrainAll(indx2, :);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMx = J0' \ dMr(indx, :);
    dMTildex = J0' \ dMTilder(indx, :);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dN_X_I);
    L = BMatrix(dMx);
    G = detJ0 / detJ * (F0' \ ansatzFunctionTestFunctionEnhancedStrain);
    HLeft = BMatrix(dMTildex);
    HRight = zeros(6, 3);
    HRight(1:3, 1) = repmat(gaussPointsInSkewCoordinates(1,k)*gaussPointsInSkewCoordinates(2,k), [3, 1]);
    HRight(1:3, 2) = repmat(gaussPointsInSkewCoordinates(2,k)*gaussPointsInSkewCoordinates(3,k), [3, 1]);
    HRight(1:3, 3) = repmat(gaussPointsInSkewCoordinates(3,k)*gaussPointsInSkewCoordinates(1,k), [3, 1]);
    H = [HLeft, (F0.' \ HRight)];

    % strain tensor
    epsilonVoigt = B * uN1 + H * alphaN1e;

    if ~computePostData
        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * epsilonVoigt' * C * epsilonVoigt * detJ * gaussWeight(k);

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

function [MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M, dMr)
% Computes the shape functions for the trial function for the enhanced part
% of the strain field
numberOfGausspoints = size(gaussPointsInSkewCoordinates, 2);

MTilde = (gaussPointsInSkewCoordinates.^2).' - M * (nodesInSkewCoordinates.^2).';

dmTilder = zeros(dimension*numberOfGausspoints, dimension);
if dimension == 2
    dmTilder(1:dimension:end, 1) = 2 * gaussPointsInSkewCoordinates(1, :);
    dmTilder(2:dimension:end, 2) = 2 * gaussPointsInSkewCoordinates(2, :);
elseif dimension == 3
    dmTilder(1:dimension:end, 1) = 2 * gaussPointsInSkewCoordinates(1, :);
    dmTilder(2:dimension:end, 2) = 2 * gaussPointsInSkewCoordinates(2, :);
    dmTilder(3:dimension:end, 3) = 2 * gaussPointsInSkewCoordinates(3, :);
else
    error('Dimension not implemented!');
end
dMTilder = dmTilder - dMr * (nodesInSkewCoordinates.^2).';
end

function ansatzFunctionTestFunctionEnhancedStrainAll = computeShapeFunctionsTestFunctionEnhancedStrain(dimension, pianSumiharaAnsatzFunctionsAll, easAnsatzFunctionsAll, X, dN_xi_k_I, gaussWeight)
% Computes the shape functions for the test function for the enhanced part
% of the strain field
numberOfGausspoints = size(dN_xi_k_I, 2);
sigma = pianSumiharaAnsatzFunctionsAll;
s = zeros(size(sigma));
for i=1:size(s, 2)
    delta = zeros(size(s, 1), 1);
    for m=1:i-1
        integralTop = 0;
        integralBottom = 0;
        for k=1:numberOfGausspoints
            dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
            [~,detJ] = computeJacobian(X,dN_xi_I,1e-10);
            index = (3 * dimension - 3) * (k - 1) + 1:(3 * dimension - 3) * k;
            integralTop = integralTop + sigma(index, i)' * s(index, m)*detJ*gaussWeight(k);
            integralBottom = integralBottom + s(index, m)' * s(index, m)*detJ*gaussWeight(k);
        end
        delta = delta + integralTop/integralBottom*s(:, m);
    end
    s(:, i) = sigma(:, i) - delta;
end

e = easAnsatzFunctionsAll;
eps = zeros(size(e));
for j=1:size(e, 2)
    delta = zeros(size(s, 1), 1);
    for m=1:j-1
        integralTop = 0;
        integralBottom = 0;
        for k=1:numberOfGausspoints
            dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
            [~,detJ] = computeJacobian(X,dN_xi_I,1e-10);
            index = (3 * dimension - 3) * (k - 1) + 1:(3 * dimension - 3) * k;
            integralTop = integralTop + e(index, j)' * s(index, m)*gaussWeight(k);
            integralBottom = integralBottom + s(index, m)' * s(index, m)*detJ*gaussWeight(k);
        end
        delta = delta + integralTop/integralBottom*s(:, m);
    end
    eps(:, j) = e(:, j) - delta;
end

ansatzFunctionTestFunctionEnhancedStrainAll = eps;
end

function ansatzFunctionTestFunctionEnhancedStrainAll = computeShapeFunctionsTestFunctionEnhancedStrain2(dimension, pianSumiharaAnsatzFunctionsAll, easAnsatzFunctionsAll, XNodes, dN_xi_k_I, gaussWeight, localCoordinates)
% Computes the shape functions for the test function for the enhanced part
% of the strain field
numberOfGausspoints = size(dN_xi_k_I, 2);

numberOfNodes = size(XNodes, 2);
[hourglassVector, ~] = computeHourglassVectorAndFunction(dimension, localCoordinates);
c = 1 / numberOfNodes * XNodes * hourglassVector;
for k=1:numberOfGausspoints
    index = (3 * dimension - 3) * (k - 1) + 1:(3 * dimension - 3) * k;
    easAnsatzFunction = easAnsatzFunctionsAll(index, :);
    xi = easAnsatzFunction(1, 1);
    eta = easAnsatzFunction(2,2);
    zeta = easAnsatzFunction(3,3);
    eps = easAnsatzFunction;
    eps(:, 10:12) = zeros(6,3);
    h10_2 = 1/3*(xi*c(1,1)-zeta*c(1,3));
    h10_3 = 1/3*(xi*c(1,1)-eta*c(1,2));
    h11_1 = 1/3*(eta*c(2,2)-zeta*c(3,2));
    h11_3 = 1/3*(eta*c(2,2)-xi*c(1,2));
    h12_1 = 1/3*(zeta*c(3,3)-eta*c(2,3));
    h12_2 = 1/3*(zeta*c(3,3)-xi*c(1,3));
    eps(1, 11) = xi*zeta - h11_1;
    eps(1, 12) = xi*eta - h12_1;
    eps(2, 10) = eta*zeta - h10_2;
    eps(2, 12) = xi*eta - h12_2;
    eps(3, 10) = eta*zeta - h10_3;
    eps(3, 11) = xi*zeta - h11_3;
    easAnsatzFunctionsAll(index, :) = eps;
end
ansatzFunctionTestFunctionEnhancedStrainAll = easAnsatzFunctionsAll;
end