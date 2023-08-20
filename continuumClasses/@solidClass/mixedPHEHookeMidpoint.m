function [rData, kData, elementEnergy, array] = mixedPHEHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDPHEHOOKEMIDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'mixed'-based finite element routine with linear elastic
% material.
% This is a port-Hamiltonian formulation with independent velocities, and
% strains.
% The routine is suitable for dynamic simulation where for the midpoint 
% rule is applied.
%
% CALL
% mixedPHHookeMidpoint(obj,setupObject,computePostData)
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
% mixedPHEHyperelasticMidpoint, mixedPHEHyperelasticDiscreteGradient
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

% element degree of freedom tables and more
edof = meshObject.edof;
dimension = obj.dimension;

% data for displacement dofs
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;

% mixed FE data
Phi = mixedFEObject.shapeFunctionObject.N_k_I;
epsN = mixedFEObject.qN(e,:)';
epsN1 = dofs.edAlphaN1;
epsN05 = 1/2 *(epsN + epsN1);

% material data
lambda = materialObject.lambda;
mu = materialObject.mu;

if (dimension == 1)
    if isfield(materialObject, 'E')
        E = materialObject.E;
    else
        E = mu * (3 * lambda + 2 * mu) / (lambda + mu);
    end
else
    error('Wrong dimension. (dim=1 only)')
end
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

%% Create residual and tangent

% initialize
elementEnergy.strainEnergy = 0;
Meps = zeros(1,1);
Kveps = zeros(2,1);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    if ~computePostData
            
         % mass and stiffness matrix
         Meps = Meps + Phi(k, :)' * E * Phi(k, :) * detJ * gaussWeight(k);
         Kveps = Kveps + B'* E * Phi(k, :) * detJ * gaussWeight(k);
            
         % stored strain energy
         elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1/2* epsN1'* E *epsN1 * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        sigma_V = E*epsN1;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJStruct,stressTensor,setupObject,dimension);
    end
end

if ~computePostData
    % residual
    Rv = Kveps * epsN05;
    Reps = Meps*(epsN1 - epsN) - Kveps' * (edN1-edN)' ;
    rData{1} = Rv;
    rData{2} = Reps;
    % tangent
    kData{1,1} = zeros(2,2);
    kData{1,2} = 1/2 * Kveps;
    kData{2,1} = -Kveps';
    kData{2,2} = Meps;

end

end