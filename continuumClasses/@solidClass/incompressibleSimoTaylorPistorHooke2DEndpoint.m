function [rData, kData, elementEnergy, array] = incompressibleSimoTaylorPistorHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
numberOfDisplacementDofs = size(globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% nodal dofs
edR = obj.qR(edof(e, :), 1:dimension)';
edN = obj.qN(edof(e, :), 1:dimension)';
edN1 = dofs.edN1;
% mixed FE data
internalQN1 = mixedFEObject.qN1;
numberOfInternalDofs = size(mixedFEObject.qN1, 2);
Gamma_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
indexTheta = 1:(numberOfInternalDofs / 2);
indexP = (numberOfInternalDofs / 2 + 1):numberOfInternalDofs;
% material data
nu = materialObject.nu;
E = materialObject.E;
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
switch lower(strtok(materialObject.name, 'Hooke'))
    case 'esz'
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end
% projection tensors voigt notation
Iv = [1, 1, 0]';
Pdev = eye(3) - 1 / 2 * (Iv * Iv');
% initialization elementEnergy
elementEnergy.strainEnergy = 0;

Kuu = zeros(numberOfDisplacementDofs, numberOfDisplacementDofs);
Kup = zeros(numberOfDisplacementDofs, numberOfInternalDofs/2);
Kut = zeros(numberOfDisplacementDofs, numberOfInternalDofs/2);
Ktt = zeros(numberOfInternalDofs/2, numberOfInternalDofs/2);
Ktp = zeros(numberOfInternalDofs/2, numberOfInternalDofs/2);
Fext = zeros(numberOfDisplacementDofs, 1);
uN1 = zeros(numberOfDisplacementDofs, 1);
uN1(:) = edN1 - edR;
thetaN1 = dofs.edAlphaN1(indexTheta)';
pN1 = dofs.edAlphaN1(indexP)';

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    %shape functions enhancement
    Gamma_I = Gamma_k_I(k, :);

    %nodal operatormatix
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    Bdev = Pdev * B;
    Bvol = dN_X_I(:)';

    if ~computePostData
        %tangent and rediduum
        Kuu = Kuu + (Bdev' * C * Bdev) * detJ * gaussWeight(k);
        Kut = Kut + (Bdev' * C * 1 / 2 * Iv * Gamma_I) * detJ * gaussWeight(k);
        Kup = Kup + (Bvol' * Gamma_I) * detJ * gaussWeight(k);
        Ktt = Ktt + (1 / 4 * Gamma_I' * Iv' * C * Iv * Gamma_I) * detJ * gaussWeight(k);
        Ktp = Ktp - (Gamma_I' * Gamma_I) * detJ * gaussWeight(k);
        %         Fext = Fext + Ni'*b * detJ*gaussWeight(k);
    else
        %stresses
        %displacement gradient and strain tensor
        gradU = (edN1 - edR) * dN_X_I';
        theta = Gamma_I * thetaN1;
        p = Gamma_I * pN1;
        epsilon = zeros(3, 3);
        epsilon(1:2, 1:2) = 1 / 2 * (gradU + gradU') + 1 / 2 * theta * eye(2);
        switch lower(strtok(obj.materialObject.name, 'Hooke'))
            case 'esz'
                epsilon(3, 3) = -nu / (1 - nu) * (trace(epsilon)); %ESZ
            case 'evz'
                epsilon(3, 3) = 0; %EVZ
            otherwise
                error('not implemented')
        end
        mu = E / (2 * (1 + nu));
        %             lambda = E*nu/((1+nu)*(1-2*nu));
        %stresses
        sigmaN1 = 2 * mu * (epsilon - 1 / 3 * trace(epsilon) * eye(3)) + p * eye(3);
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
Kpp = zeros(numberOfInternalDofs/2);
rData{1} = Kuu * uN1 + Kut * thetaN1 + Kup * pN1;
rData{2} = Kut' * uN1 + Ktt * thetaN1 + Ktp * pN1;
rData{3} = Kup' * uN1 + Ktp' * thetaN1 + zeros(numberOfInternalDofs/2) * pN1;
kData{1,1} = Kuu;
kData{1,2} = Kut;
kData{1,3} = Kup;
kData{2,1} = Kut';
kData{2,2} = Ktt;
kData{2,3} = Ktp;
kData{3,1} = Kup';
kData{3,2} = Ktp';
kData{3,3} = Kpp;
end
end