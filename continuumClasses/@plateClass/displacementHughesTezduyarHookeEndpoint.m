function [rData, kData, elementEnergy, array] = displacementHughesTezduyarHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% DISPLACEMENTHUGHESTEZDUYARHOOKEENDPOINT Element routine of class plateClass.
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
% displacementHughesTezduyarHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
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
% https://doi.org/10.1115/1.3157679
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

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

numberOfNodes = size(dN_xi_k_I, 3);
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% compute material matrices
% bending component
Eb = zeros(3, 3);
Eb(1, 1) = 1;
Eb(1, 2) = nu;
Eb(2, 1) = nu;
Eb(2, 2) = 1;
Eb(3, 3) = (1 - nu) / 2;
Eb = E * h^3 / (12 * (1 - nu^2)) * Eb;
% shear component
Es = eye(2, 2);
Es = 5 / 6 * E / (2 * (1 + nu)) * h * Es;

% initialize energy
elementEnergy.strainEnergy = 0;

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

% compute Jacobian
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);

% compute direction vectors
r11 = [ed(1, 2) - ed(1, 1); ed(2, 2) - ed(2, 1)];
r21 = [ed(1, 3) - ed(1, 2); ed(2, 3) - ed(2, 2)];
r31 = [ed(1, 4) - ed(1, 3); ed(2, 4) - ed(2, 3)];
r12 = [ed(1, 4) - ed(1, 1); ed(2, 4) - ed(2, 1)];
h(1) = norm(r11);
h(2) = norm(r21);
h(3) = norm(r31);
h(4) = norm(r12);
e_var1(:, 1) = r11 / h(1);
e_var1(:, 2) = r21 / h(2);
e_var1(:, 3) = r31 / h(3);
e_var2(:, 1) = r12 / h(4);
e_var2(:, 2) = -e_var1(:, 1);
e_var2(:, 3) = -e_var1(:, 2);
e_var2(:, 4) = -e_var1(:, 3);
e_var1(:, 4) = -e_var2(:, 1);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % B-Matrix
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);
    % shear component
    Bs = zeros(2, 3*numberOfNodes);
    for a = 1:numberOfNodes
        b = mod(a, 4) + 1;
        Ga = computeGa(a, k, N_k_I, e_var1, e_var2);
        Gb = computeGa(b, k, N_k_I, e_var1, e_var2);
        Bs(:, 3*b-2) = 1 / h(a) * Ga - 1 / h(b) * Gb;
        Bs(:, 3*b-1) = (e_var2(1, b) * Ga - e_var1(1, b) * Gb) / 2;
        Bs(:, 3*b) = (e_var2(2, b) * Ga - e_var1(2, b) * Gb) / 2;
    end

    if ~computePostData
        % Tangent
        % bending component
        Kb = Bb' * Eb * Bb * detJ * gaussWeight(k);
        % shear component
        Ks = Bs' * Es * Bs * detJ * gaussWeight(k);

        KXX = KXX + Kb + Ks;
    else
        % stress computation -> TODO
        %         epsilon = B * uN1(:);
        %         sigma_V = C * epsilon;
        %         sigma = voigtToMatrix(sigma_V, 'stress');
        %         stressTensor.Cauchy = sigma;
        %         array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
RX = KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end

function Ga = computeGa(a, k, N_k_I, e_var1, e_var2)
b = mod(a, 4) + 1;
alpha_a = e_var1(:, a).' * e_var2(:, a);
alpha_b = e_var1(:, b).' * e_var2(:, b);
Ga = 1 / (1 - alpha_a^2) * N_k_I(k, a) * (e_var1(:, a) - alpha_a * e_var2(:, a)) - 1 / (1 - alpha_b^2) * N_k_I(k, b) * (e_var2(:, b) - alpha_b * e_var1(:, b));
end