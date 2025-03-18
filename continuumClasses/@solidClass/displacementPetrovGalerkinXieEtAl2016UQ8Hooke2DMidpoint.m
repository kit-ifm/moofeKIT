function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinXieEtAl2016UQ8Hooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTPETROVGALERKINHOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is US-QUAD8, where US denotes the unsymmetry of
% the stiffness matrix.
%
% CALL
% displacementPetrovGalerkinHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% https://doi.org/10.1002/nme.836 (Rajendran, Liew: A novel unsymmetric 8-node plane element immune to mesh distortion under a quadratic displacement field)
%
% SEE ALSO
% easPetrovGalerkinHookeEndpoint
%
% CREATOR(S)
% Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire shape functions
N_k_I = shapeFunctionObject.N_k_I;

% acquire the derivatives of the shape functions with respect to the parametric coordinates
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

% aquire the edof (element degree of freedom table)
edof = meshObject.edof;

% aquire the number of dimensions of the simulation
dimension = obj.dimension;

% number of Gausspoints
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% gauss weights
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

% check input
numberOfNodes = size(N_k_I, 2);
if numberOfNodes ~= 8 || numberOfGausspoints ~= 9
    error('This displacement Petrov-Galerkin formulation is only suitable for 8-node serendipity elements with a 9 point gauss integration.')
end

% aquire material data
nu = materialObject.nu;
E = materialObject.E;

% compute material matrix
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'esz'
        % ESZ
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        % EVZ
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end

% spatial coordinates of the element
x = dofs.edN1;
xN = obj.qN(edof(e, :), 1:dimension).';
% material coordinates of the element
X = obj.qR(edof(e, :), 1:dimension)';

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);
J0 = X * dN0_xi_I';

% displacement vector
uN1 = x(:) - X(:);
uN = xN(:) - X(:);
uN05 = 1/2*(uN + uN1);

% initialize residual
Re = rData{1};

% initialize tangent
Ke = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% compute skew coordinates of nodal points / gauss points
nodalPointsInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, numberOfNodes);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPoints, X, J0, shapeFunctionObject);
nodalPointsInSkewCoordinates = computeSkewCoordinates(nodalPointsInLocalCoordinates, X, J0, shapeFunctionObject);

% compute metric shape functions
[~, dM_xiBar_k_I] = computeMetricShapeFunctions(obj, dimension, nodalPointsInSkewCoordinates, gaussPointsInSkewCoordinates);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    % compute the Jacobian determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute nodal operator matrix
    B = BMatrix(dN_X_I);

    index = dimension * (k - 1) + 1:dimension * k;
    dM_X_k_I = J0' \ dM_xiBar_k_I(index, :);
    L = BMatrix(dM_X_k_I);

    if ~computePostData
        % RESIDUAL
        Re = Re + (B' * C * L) * uN05 * detJ * gaussWeight(k);

        % TANGENT
        Ke = Ke + 1/2 * B' * C * L * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (L * uN1)' * C * (L * uN1) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * B * uN1;
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = Re;

    % pass tangent
    kData{1, 1} = Ke;
end
end