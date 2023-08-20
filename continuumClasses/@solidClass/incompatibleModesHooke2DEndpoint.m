function [rData, kData, elementEnergy, array] = incompatibleModesHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% INCOMPATIBLEMODESHOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The formulation employs incompatibleModes.
%
% CALL
% incompatibleModesHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
% easHooke2DEndpoint
%
% CREATOR(S)
% Marlon Franke, Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;
BubbledM_xi_k_I = mixedFEObject.shapeFunctionObject.dM;
modificationTaylor = strcmpi(obj.elementDisplacementType, 'incompatibleModesTaylor');

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;

dimension = obj.dimension;

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

% material data
nu = materialObject.nu;
E = materialObject.E;
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'esz'
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end

% aquire the nodal values of the variables for the current element
edN1 = dofs.edN1;
edN = obj.qN(edof(e, :), 1:dimension).';
edR = obj.qR(edof(e, :), 1:dimension).';
uN1 = edN1(:) - edR(:);
uN1eBubble = dofs.edAlphaN1.';

% Jacobian matrices
J0 = edR * dN0_xi_I';
detJ0 = det(J0);

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, J, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    % bubble modes
    BubbledM_xi_I = reshape(BubbledM_xi_k_I(:,k,:),[size(BubbledM_xi_k_I,1),size(BubbledM_xi_k_I,3)]);
    if modificationTaylor
        dMxBubble = detJ0 / detJ * (J0' \ BubbledM_xi_I);
    else
        dMxBubble = J' \ BubbledM_xi_I;
    end

    % nodal operator matrix & approximation matrix
    BCompatible = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    BBubble = BMatrix(dMxBubble, 'mapVoigtObject', mapVoigtObject);

    if ~computePostData
        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * ((BCompatible * uN1(:))' * C * (BCompatible * uN1(:)) + (BBubble * uN1eBubble(:))' * C * (BBubble * uN1eBubble(:))) * detJ * gaussWeight(k);

        % TANGENT
        KDD = KDD + (BCompatible' * C * BCompatible) * detJ * gaussWeight(k);
        KDA = KDA + (BCompatible' * C * BBubble) * detJ * gaussWeight(k);
        KAD = KAD + (BBubble' * C * BCompatible) * detJ * gaussWeight(k);
        KAA = KAA + (BBubble' * C * BBubble) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        % displacement gradient and strain tensor
        U = edN1 - edR;
        gradU = U * dN_X_I';
        gradUBubble = reshape(uN1eBubble, 2, 2) * dMxBubble';
        epsilon = zeros(3, 3);
        epsilon(1:2, 1:2) = 1 / 2 * (gradU + gradU' + gradUBubble + gradUBubble');
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
        % stress at gausspoint
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
RD = KDD * uN1 + KDA * uN1eBubble;
RA = KAD * uN1 + KAA * uN1eBubble;

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
