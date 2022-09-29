function [rData, kData, elementEnergy, array] = incompressibleHughesHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
numberOfDisplacementDofs = size(globalFullEdof,2)-size(mixedFEObject.globalEdof,2);
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
numberOfInternalDOFs = size(mixedFEObject.qN1,2);
NTilde_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
indexP = 1:numberOfInternalDOFs;
% material data
selectMapVoigt(mapVoigtObject,dimension,'symmetric');
nu = materialObject.nu;
E = materialObject.E;
mu = E/(2*(1+nu));
if nu~=0.5
    kappa = E/(3*(1-2*nu));
else
    kappa = 0;
end
switch lower(strtok(obj.materialObject.name,'Hooke'))
    case 'evz'
        Cdev = mu*[1,-1,0; -1,1,0; 0,0,1]; %EVZ
    otherwise
        error('not implemented')
end
% projection tensors voigt notation
I = eye(dimension);
Iv = [1,1,0]';
% initialization elementEnergy and storageFEObject
elementEnergy.strainEnergy = 0;

Kuu = zeros(numberOfDisplacementDofs,numberOfDisplacementDofs);
Kpu = zeros(numberOfInternalDOFs,numberOfDisplacementDofs);
Kpp = zeros(numberOfInternalDOFs,numberOfInternalDOFs);
Fext = zeros(numberOfDisplacementDofs,1);
%nodal points and displacements
uN1 = zeros(numberOfDisplacementDofs,1);
uN1(:) = edN1-edR;
pN1 = dofs.edAlphaN1(indexP)';

%loop over all gauss points
for k=1:numberOfGausspoints
    %shape functions enhancement
    NTilde_I = NTilde_k_I(k,:);
    %approximation of geometry and deriviates according to x
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    %nodal operatormatix
    B = BMatrix(dN_X_I,'mapVoigtObject',mapVoigtObject);
    if ~computePostData
        %tangent and rediduum
        Kuu = Kuu + (B'*Cdev*B)*detJ*gaussWeight(k);
        Kpp = Kpp + (NTilde_I'*NTilde_I)*detJ*gaussWeight(k);
        Kpu = Kpu + (NTilde_I'*Iv'*B)*detJ*gaussWeight(k);
        %         Fext = Fext + Ni'*b * detJ*gaussWeight(k);
    else
        gradU = (edN1-edR)*dN_X_I';
        p = NTilde_I*pN1;
        epsilon = zeros(3,3);
        epsilon(1:2,1:2) = 1/2*(gradU+gradU');
        switch lower(strtok(obj.materialObject.name,'Hooke'))
            case 'esz'
                epsilon(3,3) = -nu/(1-nu)*(trace(epsilon)); %ESZ
            case 'evz'
                epsilon(3,3) = 0; %EVZ
            otherwise
                error('not implemented')
        end
        %stresses
        sigmaN1 = 2*mu*(epsilon-1/3*trace(epsilon)*eye(3)) + p*eye(3);
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJStruct,stressTensor,setupObject,dimension);
    end
end
if nu==0.5
    %in truely incompressible case special form
    KDD = Kuu; KDA = Kpu'; KAD = Kpu; KAA = zeros(numberOfInternalDOFs);
    RD = Kuu*uN1+Kpu'*pN1-Fext;
    RA = Kpu*uN1;
else
    %compressible case
    KDD = Kuu; KDA = Kpu'; KAD = Kpu; KAA = 1/kappa*Kpp;
    RD = Kuu*uN1+Kpu'*pN1-Fext;
    RA = Kpu*uN1+1/kappa*Kpp*pN1;
end
if ~computePostData
    rData{1} = RD;
    rData{2} = RA;
    kData{1,1} = KDD;
    kData{1,2} = KDA;
    kData{2,1} = KAD;
    kData{2,2} = KAA;
end
end