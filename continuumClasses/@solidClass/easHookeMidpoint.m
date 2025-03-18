function [rData, kData, elementEnergy, array] = easHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASHOOKEMIDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The routine is suitable for static and dynamic simulations where for the
% latter the midpoint rule is used ('Midpoint').
%
% CALL
% easHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% -
%
% SEE ALSO
% easHookeEndpoint
%
% CREATOR(S)
% Robin Pfefferkorn, Felix Zaehringer

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

% material data
if dimension == 2
    emod = materialObject.E;
    nu = materialObject.nu;
    DMat = emod / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; ...
        nu, 1 - nu, 0; ...
        0, 0, (1 - 2 * nu) / 2];
else
    lambda = materialObject.lambda;
    mu = materialObject.mu;
    DMat = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
        lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
        lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
        0, 0, 0, mu, 0, 0; ...
        0, 0, 0, 0, mu, 0; ...
        0, 0, 0, 0, 0, mu];
end

% aquire the nodal values of the variables for the current element
alphaN = mixedFEObject.qN;
%
qR = obj.qR;
qN = obj.qN;

edAlphaN = alphaN(e, :);
edAlphaN1 = dofs.edAlphaN1;
if dimension == 2
    alphaNVoigt = [edAlphaN(1, 1:2:end); edAlphaN(1, 2:2:end)];
    alphaN1Voigt = [edAlphaN1(1, 1:2:end); edAlphaN1(1, 2:2:end)];
else
    alphaNVoigt = [edAlphaN(1, 1:3:end); edAlphaN(1, 2:3:end); edAlphaN(1, 3:3:end)];
    alphaN1Voigt = [edAlphaN1(1, 1:3:end); edAlphaN1(1, 2:3:end); edAlphaN1(1, 3:3:end)];
end
alphaN05Voigt = 1 / 2 * (alphaNVoigt + alphaN1Voigt);
edN1 = dofs.edN1;
edN = qN(edof(e, :), 1:dimension)';
edRef = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edRef;
uN = edN - edRef;
uN05 = 1 / 2 * (uN1 + uN);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(edRef, dN_xi_k_I);
J0 = edRef * dN0_xi_I';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% initialize residual
RD = rData{1};
RA = rData{2};

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
     % compute the Jacobian determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_k_I = computedN_X_I(dN_xi_k_I, J, k);

    % shape functions for the enhanced part of the strain field
    indx2 = (3*dimension - 3) * (k - 1) + 1:(3*dimension - 3) * k;
    Er = ErAll(indx2, :);

    % nodal operator matrix & approximation matrix
    B = BMatrix(dN_X_k_I);
    G = detJ0 / detJ * (F0' \ Er);

    % strain tensor
    epsilonVoigt = B * uN1(:) + G * alphaN1Voigt(:);

    if ~computePostData
        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * epsilonVoigt' * DMat * epsilonVoigt * detJ * gaussWeight(k);

        % RESIDUAL
        RD = RD + (B' * DMat * B * uN05(:) + B' * DMat * G * alphaN05Voigt(:)) * detJ * gaussWeight(k);
        RA = RA + (G' * DMat * B * uN05(:) + G' * DMat * G * alphaN05Voigt(:)) * detJ * gaussWeight(k);

        % TANGENT
        KDD = KDD + 1 / 2 * B' * DMat * B * detJ * gaussWeight(k);
        KDA = KDA + 1 / 2 * B' * DMat * G * detJ * gaussWeight(k);
        KAD = KAD + 1 / 2 * G' * DMat * B * detJ * gaussWeight(k);
        KAA = KAA + 1 / 2 * G' * DMat * G * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = DMat * epsilonVoigt;
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

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