function [rData, kData, elementEnergy, array] = easSimoRifaiHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASSIMORIFAIHOOKEENDPOINT Element routine of class axisymmetricSolidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The routine is suitable for static and dynamic simulations where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% easSimoRifaiHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type axisymmetricSolidClass,
%      e.g. axisymmetricSolidObject.
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
% https://doi.org/10.1002/nme.1620290802
%
% SEE ALSO
% -
%
% CREATOR(S)
% Felix Zaehringer

%% SETTING
elementType = 'A'; % Alternative implementation according to paper;

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
N0_k_I = shapeFunctionObject.N0_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

edof = meshObject.edof;

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;
C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, nu, 0; nu, 1 - nu, nu, 0; nu, nu, 1 - nu, 0; 0, 0, 0, (1 - 2 * nu) / 2];

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = edN1(:) - edR(:);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
J0 = edR * dN0_xi_I.';
detJ0 = det(J0);

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% F0 Matrix
F0([1, 2, 4], [1, 2, 4]) = F0Matrix(dimension, J0);
F0(3, 3) = 1;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);

    r = N_k_I(k, :) * edR(1, :).';
    r0 = N0_k_I * edR(1, :).';

    B = zeros(4, 8);
    B([1, 2, 4], :) = BMatrix(dN_X_I);
    B(3, 1:2:end) = 1 / r * N_k_I(k, :);

    if strcmp(elementType, 'A')
        E = zeros(4, 5);
        E(1, 1) = xi;
        E(2, 2) = eta;
        E(3, 3) = xi * eta;
        E(4, 4) = xi;
        E(4, 5) = eta;

        E = r0 / r * E;
    elseif strcmp(elementType, 'B')
        a0 = 1/4*[1, 1, 1, 1].';
        a1 = 1/4*[-1, 1, 1, -1].';
        a2 = 1/4*[-1, -1, 1, 1].';
        h = 1/4*[1, -1, 1, -1].';
        xiBar = 1/3*edR(1, :)*a1/(edR(1, :)*a0);
        etaBar = 1/3*edR(1, :)*a2/(edR(1, :)*a0);
        xiEtaBar = 1/9*edR(1, :)*h/(edR(1, :)*a0);
        E = zeros(4, 5);
        E(1, 1) = xi - xiBar;
        E(2, 2) = eta - etaBar;
        E(3, 3) = xi * eta - xiEtaBar;
        E(4, 4) = xi - xiBar;
        E(4, 5) = eta - etaBar;
    elseif strcmp(elementType, 'C')
        a0 = 1/4*[1, 1, 1, 1].';
        a1 = 1/4*[-1, 1, 1, -1].';
        a2 = 1/4*[-1, -1, 1, 1].';
        xiBar = 1/3*edR(1, :)*a1/(edR(1, :)*a0);
        etaBar = 1/3*edR(1, :)*a2/(edR(1, :)*a0);
        E = zeros(4, 5);
        E(1, 1) = xi - xiBar;
        E(2, 2) = eta - etaBar;
        E(3, 3) = xi * eta *detJ/(detJ0*r);
        E(4, 4) = xi - xiBar;
        E(4, 5) = eta - etaBar;
    else
        error('Not implemented yet!');
    end

    G = detJ0 / detJ * (F0' \ E);

    if ~computePostData
        % TANGENT
        KDD = KDD + 2 * pi * (B' * C * B) * r * detJ * gaussWeight(k);
        KDA = KDA + 2 * pi * (B' * C * G) * r * detJ * gaussWeight(k);
        KAD = KAD + 2 * pi * (G' * C * B) * r * detJ * gaussWeight(k);
        KAA = KAA + 2 * pi * (G' * C * G) * r * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * (B * uN1 + G * alphaN1e);
        stressTensor.Cauchy = [sigmaN1_v(1), sigmaN1_v(3), 0; sigmaN1_v(3), sigmaN1_v(2), 0; 0, 0, sigmaN1_v(4)];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
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