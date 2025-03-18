function [rData, kData, elementEnergy, array] = pianSumiharaSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% PIANSUMIHARASAINTVENANT2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a mixed finite element formulation covering nonlinear
% mechanical processes employing the Saint-Venant material model.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor.
% The mixed formulation employs the assumed stress method by Pian & Sumihara
% by approximating the 2nd Piola-Kirchhoff stress tensor.
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% pianSumiharaSaintVenant2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, alpha)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
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
% https://doi.org/10.1186/s40323-019-0133-z
%
% SEE ALSO
% pianSumiharaHookeEndpoint
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
Pr = mixedFEObject.shapeFunctionObject.M;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);
numberOfDisplacementDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);

dimension = obj.dimension;

I = eye(dimension);

% aquire material data
lambda = materialObject.lambda;
mu = materialObject.mu;

% material matrix (EVZ)
DMat = [lambda + 2 * mu, lambda, 0; ...
    lambda, lambda + 2 * mu, 0; ...
    0, 0, mu];
Dmatinv = DMat \ eye(size(DMat, 1));

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN1 = dofs.edN1;
edAlphaN1 = dofs.edAlphaN1.';

% initialize residual
RX = rData{1};
RA = rData{2};

% initialize tangent
KXX = kData{1, 1};
KXA = kData{1, 2};
KAX = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
J0 = edR * dN0_xi_I';

F0 = F0Matrix(dimension, J0);

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

    % stress approximation (approximation of 2nd PK)
    index = 3 * (k - 1) + 1:3 * k;
    P = F0 * Pr(index, :);
    SN1Approx_v = P * edAlphaN1;
    SN1Approx = voigtToMatrix(SN1Approx_v, 'stress');

    % approximated Green-Lagrange strain tensor
    EN1Approx_v = Dmatinv * SN1Approx_v;

    if ~computePostData
        % RESIDUAL
        RX = RX + BN1.' * SN1Approx_v * detJ * gaussWeight(k);
        RA = RA + P.' * (EN1_v - EN1Approx_v) * detJ * gaussWeight(k);

        % TANGENT
        % KXX
        % geometrical tangent
        A1 = dN_X_I.' * SN1Approx * dN_X_I * detJ * gaussWeight(k);
        temporaryKXX = zeros(numberOfDisplacementDofs);
        for g = 1:dimension
            temporaryKXX(g:dimension:numberOfDisplacementDofs, g:dimension:numberOfDisplacementDofs) = A1;
        end
        KXX = KXX + temporaryKXX;

        % KXA
        KXA = KXA + BN1.' * P * detJ * gaussWeight(k);

        % KAX
        KAX = KAX + P.' * BN1 * detJ * gaussWeight(k);

        % KAA
        KAA = KAA - P.' * Dmatinv * P * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (EN1_v)' * DMat * (EN1_v) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        PN1 = FN1 * SN1Approx;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RA;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 2} = KXA;
    kData{2, 1} = KAX;
    kData{2, 2} = KAA;
end
end
