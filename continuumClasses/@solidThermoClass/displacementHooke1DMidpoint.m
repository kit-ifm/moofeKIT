function [rData, kData, elementEnergy, array] = displacementHooke1DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTHOOKE1DMIDPOINT Element routine of class solidThermoClass.
%
% FORMULATION
% This is a displacement and temperature-based 1D finite element routine 
% covering mechanical processes employing an elastic, isotropic Hooke
% ('Hooke') model (linear geometric and linear material/stress-strain 
% relation). Alltogether, it is a nonlinear formulation. 
% The formulation is based on a standard, local, reduced energy 
% formulation.
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
% displacementHooke1DHeatEquationEndpoint
% displacementHooke1DHeatEquationMidpoint
%
% CREATOR(S)
% Marlon Franke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;
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
E = materialObject.E;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;
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
    % temperature
    thetaN05e = N_k_I(k, :) * thetaN05.';
    thetaN1e = N_k_I(k, :) * thetaN1.';
    thetaNe = N_k_I(k, :) * thetaN.';
    % linearized strain
    epsilonN05 = dN_X_I * uN05(:);
    epsilonN1 = dN_X_I * uN1(:);
    epsilonN = dN_X_I * uN(:);
    % free energy
    PsiMechanicsN = 1/2 * epsilonN * E * epsilonN;
    PsiThermoN = kappa*(thetaNe - thetaR - thetaNe*log(thetaNe/thetaR));
    PsiCoupledN = -beta*(thetaNe - thetaR)*epsilonN;
    PsiN = PsiMechanicsN + PsiThermoN + PsiCoupledN;
    PsiMechanicsN1 = 1/2 * epsilonN1 * E * epsilonN1;
%     PsiThermoN1 = -kappa/thetaR*(thetaN1e - thetaR)^2;  % not energy consistent!
    PsiThermoN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiCoupledN1 = -beta*(thetaN1e - thetaR)*epsilonN1;
    PsiN1 = PsiMechanicsN1 + PsiThermoN1 + PsiCoupledN1;
    % entropy
    etaN = kappa*log(thetaNe/thetaR) + beta*epsilonN;
    etaN1 = kappa*log(thetaN1e/thetaR) + beta*epsilonN1;
    % Cauchy stress
    sigmaN1 = E * epsilonN1 - beta*(thetaN1e - thetaR);
    % internal energy
    eN = PsiN + thetaNe*etaN;
    eN1 = PsiN1 + thetaN1e*etaN1;
    if ~computePostData
        % residual
        RX = RX + dN_X_I' *( E * epsilonN05 - beta*(thetaN05e - thetaR)) * detJ * gaussWeight(k);
        RT = RT + (N_k_I(k, :)'*thetaN05e*(etaN1-etaN)/DT+dN_X_I'*k0*dN_X_I*thetaN05')*detJ*gaussWeight(k);
        % tangent
        KXX = KXX + 0.5*dN_X_I'*E*dN_X_I * detJ * gaussWeight(k);
        KXT = KXT - 0.5*dN_X_I' * beta * N_k_I(k,:) * detJ * gaussWeight(k);
        KTX = KTX + N_k_I(k,:)' * thetaN05e* 1/DT * beta * dN_X_I * detJ * gaussWeight(k);
        KTT = KTT + (N_k_I(k,:)' * ((etaN1-etaN)/DT * 0.5* N_k_I(k,:) + thetaN05e*1/DT*kappa/thetaN1e* N_k_I(k,:)) + 0.5*dN_X_I'*k0*dN_X_I) * detJ * gaussWeight(k);
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + eN1 * detJ * gaussWeight(k);
        elementEnergy.internalEnergyDifference = elementEnergy.internalEnergyDifference + (eN1 - eN) * detJ * gaussWeight(k);
    else
        % stress computation
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
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
