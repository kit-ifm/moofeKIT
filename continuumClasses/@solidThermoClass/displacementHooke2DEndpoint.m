function [rData, kData, elementEnergy, array] = displacementHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTHOOKE2DENDPOINT Element routine of class solidThermoClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is US-ATFQ4, where US denotes the unsymmetry of
% the stiffness matrix.
% This formulation extends the purely mechanical formulation to the case of
% steady-state thermoelasticity.
%
% CALL
% displacementHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data.
% kData: cell-array of size [2, 2] for tangent data.
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

% aquire material data
nu = materialObject.nu;
E = materialObject.E;
alphaT = materialObject.alphaT;
kThermal = materialObject.k;
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

x = dofs.edN1;
X = obj.qR(edof(e, :), 1:dimension)';

TRef = obj.qR(edof(e, :), dimension+1);
TN1 = dofs.thetaN1.';

% compute Jacobian
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);

% displacement vector
uN1 = x(:) - X(:);

% initialize tangent
KDD = kData{1, 1};
KDT = kData{1, 2};
KTT = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute nodal operator matrix
    B = BMatrix(dN_X_I);

    if ~computePostData
        % TANGENT
        KDD = KDD + B' * C * B * detJ * gaussWeight(k);
        KDT = KDT - (B' * C * alphaT * [1; 1; 0] * N_k_I(k, :)) * detJ * gaussWeight(k);
        KTT = KTT + (dN_X_I.' * kThermal * eye(2) * dN_X_I) * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (1 / 2 * (B * uN1)' * C * (B * uN1) + (B * uN1)' * C * alphaT * [1; 1; 0] * N_k_I(k, :) * (TN1-TRef)) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * (B * uN1 - alphaT * [1; 1; 0] * N_k_I(k, :) * (TN1-TRef));
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
RD = KDD * uN1 + KDT * (TN1-TRef);
RT = KTT * (TN1-TRef);

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RD;
    rData{2} = RT;

    % pass tangent
    kData{1, 1} = KDD;
    kData{1, 2} = KDT;
    kData{2, 2} = KTT;
end
end