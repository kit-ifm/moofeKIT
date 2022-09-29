function [rData, kData, elementEnergy, array] = easHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASHOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The routine is suitable for static and dynamic simulations where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% easHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% Marlon Franke, Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_k_I = shapeFunctionObject.dN0_xi_k_I;
E_k_I = mixedFEObject.shapeFunctionObject.N_k_I;

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
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = edN1(:) - edR(:);

% compute Jacobian matrices
dN0_xi_I = reshape(dN0_xi_k_I(:,1,:),[size(dN0_xi_k_I,1),size(dN0_xi_k_I,3)]);
J0 = edR * dN0_xi_I';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    % shape functions for the enhanced part of the strain field
    indx2 = 3 * (k - 1) + 1:3 * k;
    Er = E_k_I(indx2, :);
%     E_I = E_k_I(k,:);

    % nodal operator matrix & approximation matrix
    B = BMatrix(dN_X_I);
    G = detJ0 / detJ * (F0' \ Er);

    if ~computePostData
        % TANGENT
        KDD = KDD + (B' * C * B) * detJ * gaussWeight(k);
        KDA = KDA + (B' * C * G) * detJ * gaussWeight(k);
        KAD = KAD + (G' * C * B) * detJ * gaussWeight(k);
        KAA = KAA + (G' * C * G) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        %displacement gradient and strain tensor
        gradU = (edN1 - edR) * dN_X_I';
        epsilonEnh_v = detJ0 / detJ * (F0' \ Er) * alphaN1e;
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
        %stresses
        mu = E / (2 * (1 + nu));
        lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        sigmaN1 = lambda * trace(epsilon) * eye(3) + 2 * mu * epsilon;
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
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