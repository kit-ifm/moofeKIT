function [rData, kData, elementEnergy, array] = displacementNewConstraintHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTNEWCONSTRAINTHOOKE2DMIDPOINT Element routine of class solidVelocityClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% In this formulation, the relationship between displacements and
% velocities is established in a special way ('NewConstraint').
% The routine is suitable for dynamic simulation where the midpoint rule is
% used ('Midpoint').
%
% CALL
% displacementNewConstraintHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidVelocityClass,
%      e.g. solidVelocityObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data of every field, here:
%        X, V
% kData: cell-array of size [2, 2] for tangent data of every field, here:
%        X, V
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
meshObject = obj.meshObject;

% aquire shape functions
N_k_I = shapeFunctionObject.N_k_I;

% acquire the derivatives of the shape functions with respect to the parametric coordinates
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% aquire the edof (element degree of freedom table)
edof = meshObject.edof;

% aquire the number of dimensions of the simulation
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
rho = materialObject.rhoNew;
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

% aquire nodal positions / displacements
edN1 = dofs.edN1;
edN = obj.qN(edof(e, :), 1:dimension)';
edR = obj.qR(edof(e, :), 1:dimension)';

uN1 = edN1 - edR;
uN = edN - edR;
uN05 = 1 / 2 * (uN1 + uN);

% aquire nodal velocities
vN1 = dofs.vN1;
vN = obj.qN(edof(e, :), dimension+1:end)';
vR = obj.qR(edof(e, :), dimension+1:end)';
vN05 = 1 / 2 * (vN1 + vN);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% initialize residual
RX = rData{1};
RV = rData{2};

% initialize tangent
KXX = kData{1, 1};
KXV = kData{1, 2};
KVX = kData{2, 1};
KVV = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;
elementEnergy.kineticEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute nodal operator matrix
    B = BMatrix(dN_X_I);
    NMat = zeros(2, 2*size(N_k_I, 2));
    NMat(1, 1:2:end) = N_k_I(k, :);
    NMat(2, 2:2:end) = N_k_I(k, :);

    if ~computePostData
        % RESIDUAL
        RX = RX + (NMat.' * rho * NMat * 1/deltaT *(vN1(:) - vN(:)) + B' * C * B * uN05(:)) * detJ * gaussWeight(k);
        RV = RV + (B.' * C * B * (1/deltaT *(uN1(:) - uN(:)) - vN05(:))) * detJ * gaussWeight(k);

        % TANGENT
        KXX = KXX + 1/2 * B' * C * B * detJ * gaussWeight(k);
        KXV = KXV + NMat.' * rho * NMat * 1/deltaT * detJ * gaussWeight(k);
        KVX = KVX + B.' * C * B * 1/deltaT * detJ * gaussWeight(k);
        KVV = KVV - 1/2 * B.' * C * B * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (1 / 2 * (B * uN1(:))' * C * (B * uN1(:))) * detJ * gaussWeight(k);
        elementEnergy.kineticEnergy = elementEnergy.kineticEnergy + (1 / 2 * (NMat * vN1(:))' * rho * (NMat * vN1(:))) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * B * uN1(:);
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RV;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 2} = KXV;
    kData{2, 1} = KVX;
    kData{2, 2} = KVV;
end
end