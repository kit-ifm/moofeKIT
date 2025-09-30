function [rData, kData, elementEnergy, array] = easPetrovGalerkinNewConstraintHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASPETROVGALERKINNEWCONSTRAINTHOOKE2DMIDPOINT Element routine of class solidVelocityClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin').
% In this formulation, the relationship between displacements and
% velocities is established in a special way ('NewConstraint').
% The routine is suitable for dynamic simulation where the midpoint rule is
% used ('Midpoint').
%
% CALL
% easPetrovGalerkinNewConstraintHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidVelocityClass,
%      e.g. solidVelocityObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, V, A, B)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, V, A, B)
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
% -
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
ErAll = mixedFEObject.shapeFunctionObject.M;

edof = meshObject.edof;

dimension = obj.dimension;

% number of Gausspoints
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% gauss points / weights
gaussWeight = shapeFunctionObject.gaussWeight;

% time step size
deltaT = setupObject.timeStepSize;

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
rho = materialObject.rhoNew;

% aquire the nodal values of the variables for the current element
X = obj.qR(edof(e, :), 1:dimension).';
x = dofs.edN1;


edN1 = dofs.edN1;
edN = obj.qN(edof(e, :), 1:dimension)';
edR = obj.qR(edof(e, :), 1:dimension)';

% displacement vector
uN1 = edN1 - edR;
uN = edN - edR;
uN05 = 1 / 2 * (uN1 + uN);

% velocity vector
vN1 = dofs.vN1;
vN = obj.qN(edof(e, :), dimension+1:end)';
vN05 = 1 / 2 * (vN1 + vN);

% additional dofs vector
alphaNe = mixedFEObject.qN(e, 1:4).';
alphaN1e = dofs.edAlphaN1(1:4).';
alphaN05e = 1 / 2 * (alphaN1e + alphaNe);

betaNe = mixedFEObject.qN(e, 5:8).';
betaN1e = dofs.edAlphaN1(5:8).';
betaN05e = 1 / 2 * (betaN1e + betaNe);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);
J0 = X * dN0_xi_I';
detJ0 = det(J0);

% compute F0 matrix
F0 = F0Matrix(dimension, J0);

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, X, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, X, J0, shapeFunctionObject);
[M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[~, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr);

% initialize residual
RX = rData{1};
RV = rData{2};
RA = rData{3};
RB = rData{4};

% initialize tangent
KXX = kData{1, 1};
KXV = kData{1, 2};
KXA = kData{1, 3};
KXB = kData{1, 4};
KVX = kData{2, 1};
KVV = kData{2, 2};
KVA = kData{2, 3};
KVB = kData{2, 4};
KAX = kData{3, 1};
KAV = kData{3, 2};
KAA = kData{3, 3};
KAB = kData{3, 4};
KBX = kData{4, 1};
KBV = kData{4, 2};
KBA = kData{4, 3};
KBB = kData{4, 4};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;
elementEnergy.kineticEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % shape functions for the enhanced part of the strain field
    indx2 = 3 * (k - 1) + 1:3 * k;
    Er = ErAll(indx2, :);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMx = J0' \ dMr(indx, :);
    dMTildex = J0' \ dMTilder(indx, :);

    % nodal operator matrix & approximation matrices
    B = BMatrix(dN_X_I);
    L = BMatrix(dMx);
    G = 1 / detJ * (F0' \ Er);
    H = BMatrix(dMTildex);

    NMat = zeros(2, 2*size(N_k_I, 2));
    NMat(1, 1:2:end) = N_k_I(k, :);
    NMat(2, 2:2:end) = N_k_I(k, :);

    MMat = zeros(2, 2*size(N_k_I, 2));
    MMat(1, 1:2:end) = M_k_I(k, :);
    MMat(2, 2:2:end) = M_k_I(k, :);

    if ~computePostData
        % RESIDUAL
        RX = RX + (NMat.' * rho * NMat * 1/deltaT *(vN1(:) - vN(:)) + B' * C * L * uN05(:) + B' * C * H * alphaN05e(:)) * detJ * gaussWeight(k);
        RV = RV + (L.' * C * (L * 1/deltaT *(uN1(:) - uN(:)) + H * 1/deltaT *(alphaN1e(:) - alphaNe(:)) - B * vN05(:) - G * betaN05e(:))) * detJ * gaussWeight(k);
        RA = RA + (G.' * C * (L * uN05(:) + H * alphaN05e(:))) * detJ * gaussWeight(k);
        RB = RB + (H.' * C * (L * 1/deltaT *(uN1(:) - uN(:)) + H * 1/deltaT *(alphaN1e(:) - alphaNe(:)) - B * vN05(:) - G * betaN05e(:))) * detJ * gaussWeight(k);

        % TANGENT
        KXX = KXX + 1/2 * (B' * C * L) * detJ * gaussWeight(k);
        KXV = KXV + 1/deltaT * (NMat.' * rho * NMat) * detJ * gaussWeight(k);
        KXA = KXA + 1/2 * (B' * C * H) * detJ * gaussWeight(k);

        KVX = KVX + 1/deltaT * L.' * C * L * detJ * gaussWeight(k);
        KVV = KVV - 1/2 * L.' * C * B * detJ * gaussWeight(k);
        KVA = KVA + 1/deltaT * L.' * C * H * detJ * gaussWeight(k);
        KVB = KVB - 1/2 * (L' * C * G) * detJ * gaussWeight(k);

        KAX = KAX + 1/2 * (G' * C * L) * detJ * gaussWeight(k);
        KAA = KAA + 1/2 * (G' * C * H) * detJ * gaussWeight(k);

        KBX = KBX + 1/deltaT * (H.' * C * L) * detJ * gaussWeight(k);
        KBV = KBV - 1/2 * (H.' * C * B) * detJ * gaussWeight(k);
        KBA = KBA + 1/deltaT * (H.' * C * H) * detJ * gaussWeight(k);
        KBB = KBB - 1/2 * (H.' * C * G) * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (1 / 2 * (L * uN1(:) + H * alphaN1e(:))' * C * (L * uN1(:) + H * alphaN1e(:))) * detJ * gaussWeight(k);
        elementEnergy.kineticEnergy = elementEnergy.kineticEnergy + (1 / 2 * (NMat * vN1(:))' * rho * (NMat * vN1(:))) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        epsilonVoigt = L * uN1 + H * alphaN1e; % strain tensor
        sigmaVoigt = C * epsilonVoigt;
        sigma = voigtToMatrix(sigmaVoigt, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RV;
    rData{3} = RA;
    rData{4} = RB;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 2} = KXV;
    kData{1, 3} = KXA;
    kData{1, 4} = KXB;
    kData{2, 1} = KVX;
    kData{2, 2} = KVV;
    kData{2, 3} = KVA;
    kData{2, 4} = KVB;
    kData{3, 1} = KAX;
    kData{3, 2} = KAV;
    kData{3, 3} = KAA;
    kData{3, 4} = KAB;
    kData{4, 1} = KBX;
    kData{4, 2} = KBV;
    kData{4, 3} = KBA;
    kData{4, 4} = KBB;
end
end

function [MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr)
% Computes the shape functions for the trial function for the enhanced part
% of the strain field
MTilde = (gaussPointsInSkewCoordinates.^2).' - M_k_I * (nodesInSkewCoordinates.^2).';

dmTilder = zeros(8, 2);
dmTilder(1:2:end, 1) = 2 * gaussPointsInSkewCoordinates(1, :);
dmTilder(2:2:end, 2) = 2 * gaussPointsInSkewCoordinates(2, :);
dMTilder = dmTilder - dMr * (nodesInSkewCoordinates.^2).';
end