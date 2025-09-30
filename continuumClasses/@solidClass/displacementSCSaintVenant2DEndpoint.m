function [rData, kData, elementEnergy, array] = displacementSCSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCSAINTVENANT2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing the Saint-Venant material model.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementSCSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% -
%
% SEE ALSO
% displacementSCSaintVenantEndpoint
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

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

globalFullEdof = meshObject.globalFullEdof;
numberOfDOFs = size(globalFullEdof, 2);

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

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN1 = dofs.edN1;

% initialize residual
RX = rData{1};

% initialize tangent
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % Deformation gradient
    FN1 = edN1 * dN_X_I.';

    % B matrix
    BN1 = BMatrix(dN_X_I, FN1);

    % Cauchy-Green strain tensor
    CN1 = FN1.' * FN1;

    % Green-Lagrange strain tensor
    EN1 = 0.5 * (CN1 - I);
    EN1_v = matrixToVoigt(EN1, 'strain');

    % 2. Piola-Kirchhoff stress tensor
    SN1_v = DMat * EN1_v;
    SN1 = voigtToMatrix(SN1_v, 'stress');

    if ~computePostData
        % RESIDUAL
        RX = RX + BN1.' * SN1_v * detJ * gaussWeight(k);

        % TANGENT
        % geometrical tangent
        A1 = dN_X_I.' * SN1 * dN_X_I * detJ * gaussWeight(k);
        temporaryKXX = zeros(numberOfDOFs);
        for g = 1:dimension
            temporaryKXX(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        KXX = KXX + temporaryKXX;
        % material tangent
        KXX = KXX + BN1.' * DMat * BN1 * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (EN1_v)' * DMat * (EN1_v) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
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
