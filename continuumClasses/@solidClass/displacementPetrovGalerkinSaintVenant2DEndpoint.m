function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTPETROVGALERKINSAINTVENANT2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing the Saint-Venant material model.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
%
% CALL
% displacementPetrovGalerkinSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [1,1] for residual data of every field, here: X
% kData: cell-array of size [1, 1] for tangent data of every field, here: X
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
% displacementSCSaintVenant2DEndpoint
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

% number of Gausspoints
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% gauss weights
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

% aquire the edof (element degree of freedom table)
edof = meshObject.edof(e, :);

% aquire the number of dimensions of the simulation
dimension = obj.dimension;

% aquire the number of degrees of freedom (DOF)
globalFullEdof = meshObject.globalFullEdof;
numberOfDOFs = size(globalFullEdof, 2);

% unit matrix
I = eye(dimension);

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% material matrix
switch lower(extractAfter(obj.materialObject.name, 'SaintVenant'))
    case 'esz'
        % ESZ
        DMat = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        % EVZ
        DMat = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end

% material coordinates of the element
edR = obj.qR(edof, 1:dimension).';

% spatial coordinates of the element
edN1 = dofs.edN1;

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
J0 = edR * dN0_xi_I';

% compute skew coordinates of nodal points / gauss points
numberOfNodes = size(N_k_I, 2);
nodalPointsInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, numberOfNodes);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPoints, edR, J0, shapeFunctionObject);
nodalPointsInSkewCoordinates = computeSkewCoordinates(nodalPointsInLocalCoordinates, edR, J0, shapeFunctionObject);

% compute metric shape functions
[~, dM_xiBar_k_I] = computeMetricShapeFunctions(obj, dimension, nodalPointsInSkewCoordinates, gaussPointsInSkewCoordinates);

% initialize residual
RX = rData{1};

% initialize tangent
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    
    index = dimension * (k - 1) + 1:dimension * k;
    dM_X_I = J0.' \ dM_xiBar_k_I(index, :);

    % Deformation gradient
    FHatN1 = edN1 * dM_X_I.';

    % Nodal operator matrices
    BN1 = BMatrix(dN_X_I, FHatN1);
    LN1 = BMatrix(dM_X_I, FHatN1);

    % Cauchy-Green strain tensor
    CHatN1 = FHatN1.' * FHatN1;

    % Green-Lagrange strain tensor
    EHatN1 = 0.5 * (CHatN1 - I);
    EHatN1_v = matrixToVoigt(EHatN1, 'strain');

    % 2. Piola-Kirchhoff stress tensor
    SHatN1_v = DMat * EHatN1_v;
    SHatN1 = voigtToMatrix(SHatN1_v, 'stress');

    if ~computePostData
        % RESIDUAL
        RX = RX + BN1.' * SHatN1_v * detJ * gaussWeight(k);

        % TANGENT
        % geometrical tangent
        A1 = dN_X_I.' * SHatN1 * dM_X_I * detJ * gaussWeight(k);
        temporaryKXX = zeros(numberOfDOFs);
        for g = 1:dimension
            temporaryKXX(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        KXX = KXX + temporaryKXX;
        % material tangent
        KXX = KXX + BN1.' * DMat * LN1 * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (EHatN1_v)' * DMat * (EHatN1_v) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        PHatN1 = FHatN1 * SHatN1;
        stressTensor.FirstPK = PHatN1;
        stressTensor.Cauchy = 1 / det(FHatN1) * PHatN1 * FHatN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;

    % pass tangent
    kData{1, 1} = KXX;
end
end
