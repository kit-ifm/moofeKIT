function [rData, kData, elementEnergy, array] = easPetrovGalerkinHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASPETROVGALERKINHOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is Q1U/E4, where U denotes the unsymmetry of
% the stiffness matrix.
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
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, C, G, c, LambdaC, LambdaG, Lambdac)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, C, G, c, LambdaC,
%        LambdaG, Lambdac)
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
% easHooke2DEndpoint
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
J = X * dNr';
JN1 = x * dNr';
J0 = X * dNr0';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, X, J0);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, X, J0);
[M, dMr] = computeShapeFunctionsTrialFunctionDisplacement(nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[~, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M, dMr);

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
    indx2 = 3 * (k - 1) + 1:3 * k;
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

    if ~computePostData
        % TANGENT
        KDD = KDD + (B' * C * L) * detJ * gaussWeight(k);
        KDA = KDA + (B' * C * H) * detJ * gaussWeight(k);
        KAD = KAD + (G' * C * L) * detJ * gaussWeight(k);
        KAA = KAA + (G' * C * H) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        % displacement gradient and strain tensor
        gradU = (x - X) * dNx';
        epsilonEnh_v = G * alphaN1e;
        epsilon = zeros(3, 3);
        epsilon(1:2, 1:2) = 1 / 2 * (gradU + gradU') + epsilonEnh_v([1, 3; 3, 2]);
        switch lower(strtok(obj.materialObject.name, 'Hooke'))
            case 'esz'
                epsilon(3, 3) = -nu / (1 - nu) * (trace(epsilon)); %ESZ
            case 'evz'
                epsilon(3, 3) = 0; %EVZ
            otherwise
                error('not implemented')
        end
        % stresses
        mu = E / (2 * (1 + nu));
        lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        sigmaN1 = lambda * trace(epsilon) * eye(3) + 2 * mu * epsilon;
        stressTensor.Cauchy = sigmaN1;
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

function xiBar = computeSkewCoordinates(xi, XNodes, J0)
% Computes the skew coordinates for given local coordinates
H1 = xi(1, :) .* xi(2, :);
h1 = [1; -1; 1; -1];
c1 = 1 / 4 * XNodes * h1;
xiBar = xi + (J0 \ c1) * H1;
end

function [M, dMr] = computeShapeFunctionsTrialFunctionDisplacement(nodesInSkewCoordinates, gaussPointsInSkewCoordinates)
% Computes the shape functions for the trial function for the displacement
numberOfGausspoints = size(gaussPointsInSkewCoordinates, 2);
A = [ones(4, 1), nodesInSkewCoordinates(1, :).', nodesInSkewCoordinates(2, :).', (nodesInSkewCoordinates(1, :) .* nodesInSkewCoordinates(2, :)).'];
B = [ones(numberOfGausspoints, 1), gaussPointsInSkewCoordinates(1, :).', gaussPointsInSkewCoordinates(2, :).', (gaussPointsInSkewCoordinates(1, :) .* gaussPointsInSkewCoordinates(2, :)).'];
M = B / A;

C = zeros(2*numberOfGausspoints, 4);
C(1:2:end, 2) = 1;
C(2:2:end, 3) = 1;
C(1:2:end, 4) = gaussPointsInSkewCoordinates(2, :).';
C(2:2:end, 4) = gaussPointsInSkewCoordinates(1, :).';
dMr = C / A;
end

function [MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M, dMr)
% Computes the shape functions for the trial function for the enhanced part
% of the strain field
MTilde = (gaussPointsInSkewCoordinates.^2).' - M * (nodesInSkewCoordinates.^2).';

dmTilder = zeros(8, 2);
dmTilder(1:2:end, 1) = 2 * gaussPointsInSkewCoordinates(1, :);
dmTilder(2:2:end, 2) = 2 * gaussPointsInSkewCoordinates(2, :);
dMTilder = dmTilder - dMr * (nodesInSkewCoordinates.^2).';
end