function [rData, kData, elementEnergy, array] = pianSumiharaHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% Pian-Sumihara element ansatz

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfDisplacementDofs = size(globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;
% nodal dofs
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
% mixed FE data
Pr = mixedFEObject.shapeFunctionObject.M;
numberOfInternalDofs = size(mixedFEObject.qN1, 2);
% material data and voigt notation
nu = materialObject.nu;
E = materialObject.E;
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
switch lower(strtok(materialObject.name, 'Hooke'))
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2]; %EVZ
    otherwise
        error('not implemented')
end
Cinv = C \ eye(size(C, 1));
% initialize elementEnergy
elementEnergy.strainEnergy = 0;

Ge = zeros(numberOfInternalDofs, numberOfDisplacementDofs);
He = zeros(numberOfInternalDofs, numberOfInternalDofs);
edN1 = dofs.edN1;
edN = qN(edof(e, :), 1:dimension)';
edR = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edR;
betaN1e = dofs.edAlphaN1';

%jacobian element centroid
J0 = edR * dN0_xi_I';
F0 = F0Matrix(dimension, J0);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    %nodal operatormatix
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);

    %stress approximation
    index = 3 * (k - 1) + 1:3 * k;
    P = F0 * Pr(index, :);
    if ~computePostData
        %tangent and rediduum
        Ge = Ge + P' * B * detJ * gaussWeight(k);
        He = He - P' * Cinv * P * detJ * gaussWeight(k);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1(:))' * C * (B * uN1(:)) * detJ * gaussWeight(k);
    else
        %stresses
        sigmaPSN1 = F0 * Pr(index, :) * betaN1e;
        sigmaN1 = zeros(3, 3);
        sigmaN1(1:2, 1:2) = sigmaPSN1([1, 3; 3, 2]);
        sigmaN1(3, 3) = nu * trace(sigmaN1);
        % stress at gausspoint
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = Ge' * betaN1e;
    rData{2} = Ge * uN1(:) + He * betaN1e;
    kData{1, 1} = zeros(numberOfDisplacementDofs);
    kData{1, 2} = Ge';
    kData{2, 1} = Ge;
    kData{2, 2} = He;
end
end
