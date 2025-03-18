function [rData, kData, elementEnergy, array] = incompressibleHerrmannHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = linearEndPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
% homogenous linear-elastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
%
% 18.01.2016 M.FRANKE

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
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% nodal dofs
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
% mixed FE data
numberOfInternalDofs = size(mixedFEObject.qN1, 2);
numberOfDisplacementDofs = size(globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
NTilde_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
indexP = 1:numberOfInternalDofs;
% material data
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
nu = materialObject.nu;
E = materialObject.E;
mu = E / (2 * (1 + nu));
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'evz'
        Cmu = mu * diag([2, 2, 1]); %EVZ
    otherwise
        error('not implemented')
end
% projection tensors voigt notation
Iv = [1, 1, 0]';
I = eye(3);
% initialization elementEnergy
elementEnergy.strainEnergy = 0;

Kuu = zeros(numberOfDisplacementDofs, numberOfDisplacementDofs);
Kpu = zeros(numberOfInternalDofs, numberOfDisplacementDofs);
Kpp = zeros(numberOfInternalDofs, numberOfInternalDofs);
Fext = zeros(numberOfDisplacementDofs, 1);
%nodal points and displacements
uN1 = zeros(numberOfDisplacementDofs, 1);
uN1(:) = edN1 - edR;
pN1 = dofs.edAlphaN1(indexP)';

%loop over all gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    %shape functions enhancement
    Ntilde = NTilde_k_I(k, :);
    %nodal operatormatix
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    if ~computePostData
        %tangent and rediduum
        Kuu = Kuu + (B' * Cmu * B) * detJ * gaussWeight(k);
        Kpp = Kpp + (Ntilde' * Ntilde) * detJ * gaussWeight(k);
        Kpu = Kpu + (Ntilde' * Iv' * B) * detJ * gaussWeight(k);
        %         Fext = Fext + Ni'*b * detJ*gaussWeight(k);
    else
        gradU = (edN1 - edR) * dN_X_I';
        p = 2 / 3 * mu * trace(gradU);
        epsilon = zeros(3, 3);
        epsilon(1:2, 1:2) = 1 / 2 * (gradU + gradU');
        switch lower(strtok(obj.materialObject.name, 'Hooke'))
            case 'esz'
                epsilon(3, 3) = -nu / (1 - nu) * (trace(epsilon)); %ESZ
            case 'evz'
                epsilon(3, 3) = 0; %EVZ
            otherwise
                error('not implemented')
        end
        %stresses
        sigmaN1 = 2 * mu * (epsilon - 1 / 3 * trace(epsilon) * eye(3)) + p * eye(3);
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = Kuu * uN1 + Kpu' * pN1 - Fext;
    rData{2} = Kpu * uN1 + 1 / lambda * Kpp * pN1;
    kData{1, 1} = Kuu;
    kData{1, 2} = Kpu';
    kData{2, 1} = Kpu;
    kData{2, 2} = 1 / lambda * Kpp;
end
end
