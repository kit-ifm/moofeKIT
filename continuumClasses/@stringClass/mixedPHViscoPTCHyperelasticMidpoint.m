function [rData, kData, elementData, array] = mixedPHViscoPTCHyperelasticMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDPHHYPERELASTICDISCRETEGRADIENT Element routine of class solidClass.
%
% FORMULATION
% This is a 'mixed'-based finite element routine with linear elastic
% material (St.-Venant Kirchhoff material for large displacements).
% This is a port-Hamiltonian formulation with independent velocities, and
% strains.
% The routine is suitable for dynamic simulation where the midpoint 
% rule is applied.
%
% CALL
% mixedPHStVenantMidpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (...)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (...)
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
% ...
%
% SEE ALSO
% mixedPHCStVenantMidpoint, mixedPHCHyperelasticMidpoint
%
% CREATOR(S)
% Philipp Kinon

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
mapVoigtObject = obj.mapVoigtObject;
meshObject = obj.meshObject;
h = setupObject.timeStepSize;
% element degree of freedom tables and more
edof = meshObject.edof;
dimension = obj.dimension;

% data for displacement dofs
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% aquire the nodal values of the variables for the current element
dimAbs = size(obj.qR,2); 
edR = obj.qR(edof(e, :), 1:dimAbs).';
edN = obj.qN(edof(e, :), 1:dimAbs).';
edN1 = dofs.edN1;

% position vectors
rN = edN(:);
rN1 = edN1(:);
rN05 = 1/2*(rN+rN1);

% mixed FE data
Phi = mixedFEObject.shapeFunctionObject.N_k_I;
CN = mixedFEObject.qN(e,1)';
CN1 = dofs.edAlphaN1(1);
CN05 = 1/2 *(CN + CN1);
CelN = mixedFEObject.qN(e,2)';
CelN1 = dofs.edAlphaN1(2);
CelN05 = 1/2 *(CelN + CelN1);

% material data
EA = materialObject.EA/2;
EAel = materialObject.EA/2;
etaA = materialObject.etaA;

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

%% Create residual and tangent
% initialize
elementData.strainEnergy = 0;
elementData.dissipatedEnergy = 0;
Ms = zeros(1,1);
R_eta = zeros(1,1);
K = zeros(2*dimAbs,1);
DHDC = 0;
DHDCel = 0;
DKDr = zeros(2*dimAbs,2*dimAbs);
D2HDC2 = 0;
D2HDCel2 = 0;
% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
% Run through all Gauss points
for k = 1:numberOfGausspoints

    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    B = [dN_X_I(1)*eye(dimAbs) dN_X_I(2)*eye(dimAbs)];

    if ~computePostData

         % mass and stiffness matrix
         Ms = Ms + Phi(k, :)' * Phi(k, :) * detJ * gaussWeight(k);
         R_eta = R_eta + 2/etaA * Phi(k, :)' * Phi(k, :) * detJ * gaussWeight(k);
         K = K + B'*B*rN05*Phi(k, :) * detJ * gaussWeight(k);
         DKDr = DKDr + B'*B*Phi(k, :) * detJ * gaussWeight(k);

         %current strains
         CN = Phi(k,:)*CN;
         CN1 = Phi(k,:)*CN1;
         CN05 = Phi(k,:)*CN05;

         % use midpoint rule 
         DHDC = DHDC + Phi(k,:)'*getStrainEnergyDerivative(EA,CN05) * detJ * gaussWeight(k);
         DHDCel = DHDCel + Phi(k,:)'*getStrainEnergyDerivative(EAel,CelN05) * detJ * gaussWeight(k);
         D2HDC2 = D2HDC2 + Phi(k,:)'*getStrainEnergyHessian(EA,CN05) * detJ * gaussWeight(k);
         D2HDCel2 = D2HDCel2 + Phi(k,:)'*getStrainEnergyHessian(EAel,CelN05) * detJ * gaussWeight(k);
         % stored strain energy
         elementData.strainEnergy = elementData.strainEnergy + getStrainEnergy(EA,CN1) * detJ * gaussWeight(k) + getStrainEnergy(EA,CelN1) * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        [~, detJStruct, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
        SN1_V = 2*getStrainEnergyDerivative(EA,CN1);
        SN1 = voigtToMatrix(SN1_V, 'stress');
        stressTensor.Cauchy = SN1;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJStruct,stressTensor,setupObject,dimension);
    end
    
    
end

MsInv = Ms \ eye(size(Ms));
elementData.dissipatedEnergy = elementData.dissipatedEnergy + h*(MsInv*DHDCel)'*R_eta*MsInv*DHDCel;

if ~computePostData
    
    % residual
    Rv = 2 * K * MsInv * DHDC + 2 * K * MsInv * DHDCel;
    Reps = Ms * (CN1 - CN) - 2 * K' * (rN1-rN) ;
%    Revo = Ms * (CelN1 - CelN) - 2 * K' * (rN1-rN) +  h * R_eta * MsInv * DHDCel; 
    Revo = Ms * (CelN1 - CelN) - 2 * K' * (rN1-rN) +  h * 2/etaA * DHDCel; 

    rData{1} = Rv;
    rData{2} = Reps;
    rData{3} = Revo;
    
    % tangent
    kData{1,1} = 2 * DKDr * MsInv * DHDC * 1/2 + 2 * DKDr * MsInv * DHDCel * 1/2;
    kData{1,2} = 2 * K * MsInv * D2HDC2 *1/2;
    kData{1,3} = 2 * K * MsInv * D2HDCel2*1/2;
    kData{2,1} = -2 * K' - 2*(DKDr*(rN1-rN))'*1/2;
    kData{2,2} = Ms;
    kData{2,3} = 0;
    kData{3,1} = -2 * K' - 2*(DKDr*(rN1-rN))'*1/2;
    kData{3,2} = 0;
    kData{3,3} = Ms + h*2/etaA*D2HDCel2 * 1/2;
    
end

end

%% Constitutive model
function W = getStrainEnergy(EA, C)

    W = 1/4*EA*(C-log(C)-1);

end

function DWC = getStrainEnergyDerivative(EA, C)

    DWC = 1/4*EA*(1-1/C);

end

function DWC = getStrainEnergyHessian(EA, C)

    DWC = 1/4*EA*1/(C^2);

end