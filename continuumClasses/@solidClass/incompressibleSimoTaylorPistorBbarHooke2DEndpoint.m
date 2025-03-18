function [rData, kData, elementEnergy, array] = incompressibleSimoTaylorPistorBbarHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
Gamma_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
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
I = eye(dimension);
Iv = [1, 1, 0]';
% initialization tangent & elementEnergy
KXX = kData{1, 1};
elementEnergy.strainEnergy = 0;

%nodal points and displacements
uN1 = zeros(numberOfDisplacementDofs, 1);
uN1(:) = edN1 - edR;

%loop over all gauss points (first loop)
He = zeros(size(Gamma_k_I, 2));
bGamma = zeros(size(Gamma_k_I, 2), size(dN_xi_k_I, 2)*dimension);
for k = 1:numberOfGausspoints
    % compute the Jacobian determinant
    [detJ, ~, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    %shape functions
    Gamma_I = Gamma_k_I(k, :);

    bGamma = bGamma + Gamma_I' * (dN_X_I(:))' * detJ * gaussWeight(k);
    He = He + Gamma_I * Gamma_I' * detJ * gaussWeight(k);
end
%loop over all gauss points (second loop)
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    %shape functions enhancement
    Gamma_I = Gamma_k_I(k, :);

    %approximation of geometry and deriviates according to x

    %Bmatrices
    bbar = bGamma * (He \ Gamma_I);
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    Bvol = 1 / 2 * Iv * (dN_X_I(:))';
    Bvolbar = 1 / 2 * Iv * (bbar(:))';
    Bbar = B - Bvol + Bvolbar;

    if ~computePostData
        %tangent and rediduum
        KXX = KXX + (Bbar' * C * Bbar) * detJ * gaussWeight(k);
        %         Fext = Fext + Ni'*b * detJ*gaussWeight(k);
    else
        %stresses
        sigmaN1v = C * Bbar * uN1;
        sigmaN1 = zeros(3, 3);
        sigmaN1(1:2, 1:2) = sigmaN1v([1, 3; 3, 2]);
        sigmaN1(3, 3) = nu * trace(sigmaN1);
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = KXX * uN1;
    kData{1, 1} = KXX;
end
end
