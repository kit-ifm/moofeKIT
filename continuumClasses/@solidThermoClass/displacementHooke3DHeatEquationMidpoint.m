function [rData, kData, elementEnergy, array] = displacementHooke3DHeatEquationMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% DISPLACEMENTHOOKE3DHEATEQUATIONMIDPOINT Element routine of class 
% solidThermoClass.
%
% FORMULATION
% This is a displacement and temperature-based 1D finite element routine 
% covering mechanical processes employing an elastic, isotropic Hooke
% ('Hooke') model (linear geometric and linear material/stress-strain 
% relation). Alltogether, it is a nonlinear formulation. 
% The formulation is based on a heat equation instead of standard, local, 
% reduced energy formulation.
% Implementation is due to work-conjugated Cauchy stress tensor and 
% linearized strain tensor.
% The routine is suitable for dynamic simulation due to midpoint
% ('Midpoint') integration scheme and should be energy consistent for 
% suitable free energy formulations.
%
% CALL
% displacementHooke1DMidpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, T)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, T)
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
% https://doi.org/10.1007/978-90-481-2331-5
% 
% SEE ALSO  
% displacementHooke1DEndpoint
% displacementHooke1DMidpoint
% displacementHooke3DHeatEquationEndpoint
%
% CREATOR(S)
% Marlon Franke, Lunes Hamadi

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;
mapVoigtObject = obj.mapVoigtObject;
% shape functions
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% 
dimension = obj.dimension;
edof = meshObject.edof;
DT = setupObject.timeStepSize;
% material data
kappa = materialObject.kappa;
alphaTensor = materialObject.alpha;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;
rho = materialObject.rho;
lambda = materialObject.lambda;
mu = materialObject.mu;

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
if (dimension == 1)
    if isfield(materialObject, 'E')
        C = materialObject.E;
    else
        C = mu * (3 * lambda + 2 * mu) / (lambda + mu);
    end
else
    C = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
        lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
        lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
        0, 0, 0, mu, 0, 0; ...
        0, 0, 0, 0, mu, 0; ...
        0, 0, 0, 0, 0, mu];
end

alpha = zeros(6,1);
for i =1:3
    alpha(i,1) = alphaTensor(i,i);
end
beta = C * alpha;

% initialize sub-residual vector
RX = rData{1};
RT = rData{2};

% initialize sub-tangent matrices
KXX = kData{1, 1};
KXT = kData{1, 2};
KTX = kData{2, 1};
KTT = kData{2, 2};

% degrees of freedom
qR = obj.qR;
qN = obj.qN;
edN1 = dofs.edN1;
edN = qN(edof(e,:),1:dimension).';
edR = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edR;
uN = edN - edR;
uN05 = 0.5*(uN + uN1);
thetaN1 = dofs.thetaN1;
thetaN = qN(edof(e,:),dimension+1).';
thetaN05 = 0.5*(thetaN + thetaN1);

% initialize energy
elementEnergy.internalEnergy = 0;
elementEnergy.internalEnergyDifference = 0;

%% Gauss loop
for k = 1:numberOfGausspoints
     [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
   
     % temperature
    thetaN05GP = N_k_I(k, :) * thetaN05(:);
    thetaN1GP = N_k_I(k, :) * thetaN1(:);
    thetaNGP = N_k_I(k, :) * thetaN(:);
    thetaRGP = N_k_I(k, :) * thetaR(:);
    
    % liN_k_Iearized strain
    epsilonN1 = B * uN1(:);
    epsilonN = B * uN(:) ;

    % Cauchy stress
    sigmaN1 = C * epsilonN1 - beta*(thetaN1GP  - thetaRGP );
    % internal energy
    eN = 1 / 2 * epsilonN' * C * epsilonN + thetaRGP *beta'*epsilonN + rho*kappa*(thetaNGP -thetaRGP );
    eN1 = 1 / 2 * epsilonN1'*  C  * epsilonN1+ thetaRGP *beta'*epsilonN1 + rho*kappa*(thetaN1GP -thetaRGP );
    
    if ~computePostData
        % residual
        RX = RX + (B' * C * B* uN05(:) - B' * beta * (thetaN05GP  - thetaRGP ) ) * detJ * gaussWeight(k);
        RT = RT + (N_k_I(k, :)' * rho * kappa * (thetaN1GP  - thetaNGP ) / DT + N_k_I(k,:)'* thetaN05GP  * beta'* (epsilonN1-epsilonN)/DT + (dN_X_I' * rho * k0 * dN_X_I)*thetaN05(:)) * detJ * gaussWeight(k);
       
        % tangent
        KXX = KXX + 0.5 * B' * C * B * detJ * gaussWeight(k);
        KXT = KXT - 0.5 * B' * beta * N_k_I(k, :) * detJ * gaussWeight(k);
        KTX = KTX + N_k_I(k,:)' * beta' * thetaN05GP  / DT * B * detJ * gaussWeight(k);
        KTT = KTT + (N_k_I(k,:)' * rho* kappa / DT*N_k_I(k,:) + 0.5 * N_k_I(k,:)' * beta' * (epsilonN1 - epsilonN) / DT * N_k_I(k,:) + 0.5 * dN_X_I' * rho * k0 * dN_X_I) * detJ * gaussWeight(k);
        
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + eN1 * detJ * gaussWeight(k);
        elementEnergy.internalEnergyDifference = elementEnergy.internalEnergyDifference + (eN1 - eN) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        sigma = voigtToMatrix(sigmaN1, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N, k, gaussWeight, detStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    rData{2} = RT;
    kData{1, 1} = KXX;
    kData{1, 2} = KXT;
    kData{2, 1} = KTX;
    kData{2, 2} = KTT;
end
end