function [rData, kData, elementEnergy, array] = easPetrovGalerkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASPETROVGALERKINHOOKEENDPOINT Element routine of class solidClass.
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
% easPetrovGalerkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% easPetrovGalerkinHooke2DEndpoint
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
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
dNr0 = shapeFunctionObject.dNr0;
ErAll = mixedFEObject.shapeFunctionObject.N;

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
x = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = x(:) - X(:);

% compute Jacobian matrices
J = X * dNr';
JN1 = x * dNr';
J0 = X * dNr0';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 2^dimension);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(dimension, nodesInLocalCoordinates, X, J0);
gaussPointsInSkewCoordinates = computeSkewCoordinates(dimension, gaussPointsInLocalCoordinates, X, J0);
[M, dMr] = computeShapeFunctionsTrialFunctionDisplacement(dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[~, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M, dMr);

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    detJN1 = det(JN1(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end

    % shape functions for the enhanced part of the strain field
    indx2 = (3 * dimension - 3) * (k - 1) + 1:(3 * dimension - 3) * k;
    Er = ErAll(indx2, :);

    % derivatives with respect to the physical coordinates
    dNx = J(:, indx)' \ dNr(indx, :);
    dMx = J0' \ dMr(indx, :);
    dMTildex = J0' \ dMTilder(indx, :);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dNx);
    L = BMatrix(dMx);
    G = detJ0 / detJ * (F0' \ Er);
    H = BMatrix(dMTildex);

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
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
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

function xiBar = computeSkewCoordinates(dimension, localCoordinates, XNodes, J0)
% Computes the skew coordinates for given local coordinates
numberOfNodes = size(XNodes, 2);
[hourglassVector, hourglassFunction] = computeHourglassVectorAndFunction(dimension, localCoordinates);

H1 = localCoordinates(1, :) .* localCoordinates(2, :);
c = 1 / numberOfNodes * XNodes * hourglassVector;
xiBar = localCoordinates + (J0 \ c) * H1;
end

function [hourglassVector, hourglassFunction] = computeHourglassVectorAndFunction(dimension, localCoordinates)
xi = localCoordinates(1, :);
eta = localCoordinates(2, :);
if dimension == 2
    hourglassVector = [+1; -1; +1; -1];
    hourglassFunction = xi .* eta;
elseif dimension == 3
    hourglassVector = zeros(8, 4);
    hourglassVector(:, 1) = [+1; +1; -1; -1; -1; -1; +1; +1];
    hourglassVector(:, 2) = [+1; -1; -1; +1; -1; +1; +1; -1];
    hourglassVector(:, 3) = [+1; -1; +1; -1; +1; -1; +1; -1];
    hourglassVector(:, 4) = [-1; +1; -1; +1; +1; -1; +1; -1];

    zeta = localCoordinates(3, :);
    hourglassFunction = zeros(size(localCoordinates, 2), 4);
    hourglassFunction(1, :) = eta .* zeta;
    hourglassFunction(2, :) = xi .* zeta;
    hourglassFunction(3, :) = xi .* eta;
    hourglassFunction(4, :) = xi .* eta .* zeta;
else
    error('Dimension not implemented!');
end
end

function [M, dMr] = computeShapeFunctionsTrialFunctionDisplacement(dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates)
% Computes the shape functions for the trial function for the displacement
numberOfNodes = size(nodesInSkewCoordinates, 2);
numberOfGausspoints = size(gaussPointsInSkewCoordinates, 2);
xiBarNodes = nodesInSkewCoordinates(1, :).';
etaBarNodes = nodesInSkewCoordinates(2, :).';
xiBarGaussPoints = gaussPointsInSkewCoordinates(1, :).';
etaBarGaussPoints = gaussPointsInSkewCoordinates(2, :).';
if dimension == 2
    A = [ones(numberOfNodes, 1), xiBarNodes, etaBarNodes, xiBarNodes .* etaBarNodes];
    B = [ones(numberOfGausspoints, 1), xiBarGaussPoints, etaBarGaussPoints, xiBarGaussPoints .* etaBarGaussPoints];
elseif dimension == 3
    zetaBarNodes = nodesInSkewCoordinates(3, :).';
    zetaBarGaussPoints = gaussPointsInSkewCoordinates(3, :).';
    A = [ones(numberOfNodes, 1), xiBarNodes, etaBarNodes, zetaBarNodes, etaBarNodes .* zetaBarNodes, zetaBarNodes .* xiBarNodes, xiBarNodes .* etaBarNodes, xiBarNodes .* etaBarNodes .* zetaBarNodes];
    B = [ones(numberOfGausspoints, 1), xiBarGaussPoints, etaBarGaussPoints, zetaBarGaussPoints, etaBarGaussPoints .* zetaBarGaussPoints, zetaBarGaussPoints .* xiBarGaussPoints, xiBarGaussPoints .* etaBarGaussPoints, xiBarGaussPoints .* etaBarGaussPoints .* zetaBarGaussPoints];
else
    error('Dimension not implemented!');
end
M = B / A;

C = zeros(dimension*numberOfGausspoints, numberOfNodes);
if dimension == 2
    C(1:dimension:end, 2) = 1;
    C(2:dimension:end, 3) = 1;
    C(1:dimension:end, 4) = etaBarGaussPoints;
    C(2:dimension:end, 4) = xiBarGaussPoints;
elseif dimension == 3
    C(1:dimension:end, 2) = 1;
    C(2:dimension:end, 3) = 1;
    C(3:dimension:end, 4) = 1;
    C(2:dimension:end, 5) = zetaBarGaussPoints;
    C(3:dimension:end, 5) = etaBarGaussPoints;
    C(1:dimension:end, 6) = zetaBarGaussPoints;
    C(3:dimension:end, 6) = xiBarGaussPoints;
    C(1:dimension:end, 7) = etaBarGaussPoints;
    C(2:dimension:end, 7) = xiBarGaussPoints;
    C(1:dimension:end, 8) = etaBarGaussPoints .* zetaBarGaussPoints;
    C(2:dimension:end, 8) = xiBarGaussPoints .* zetaBarGaussPoints;
    C(3:dimension:end, 8) = xiBarGaussPoints .* etaBarGaussPoints;
else
    error('Dimension not implemented!');
end
dMr = C / A;
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