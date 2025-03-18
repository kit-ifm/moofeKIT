function [rData, kData, elementEnergy, array] = displacementHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTHOOKEENDPOINT Element routine of class axisymmetricSolidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The routine is suitable for static and dynamic simulations where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type axisymmetricSolidClass,
%      e.g. axisymmetricSolidObject.
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

% gauss weights
gaussWeight = shapeFunctionObject.gaussWeight;

% aquire material data
nu = materialObject.nu;
E = materialObject.E;

% compute material matrix
C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0, nu; nu, 1 - nu, 0, nu; 0, 0, (1 - 2 * nu) / 2, 0; nu, nu, 0, 1 - nu];

% spatial coordinates of the element
x = dofs.edN1;
% material coordinates of the element
X = obj.qR(edof(e, :), 1:dimension)';

% compute Jacobian
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);

% displacement vector
uN1 = x(:) - X(:);

% initialize residual
Re = rData{1};

% initialize tangent
Ke = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    % compute the Jacobian determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    r = N_k_I(k, :) * X(1, :).';

    % compute nodal operator matrix
    B = zeros(4, length(uN1));
    B(1:3, :) = BMatrix(dN_X_I);
    B(4, 1:2:end) = 1 / r * N_k_I(k, :);

    if ~computePostData
        % RESIDUAL
        Re = Re + 2 * pi * (B' * C * B) * r * uN1 * detJ * gaussWeight(k);

        % TANGENT
        % Note: Corresponds to the elements stiffness matrix
        Ke = Ke + 2 * pi * B' * C * B * r * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1)' * C * (B * uN1) * r * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * B * uN1;
        stressTensor.Cauchy = [sigmaN1_v(1), sigmaN1_v(3), 0; sigmaN1_v(3), sigmaN1_v(2), 0; 0, 0, sigmaN1_v(4)];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
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