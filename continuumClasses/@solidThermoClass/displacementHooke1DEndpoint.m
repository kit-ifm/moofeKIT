function [rData, kData, elementEnergy, array] = displacementHooke1DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTHOOKE1DENDPOINT Element routine of class solidThermoClass.
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
% The routine is suitable for static and dynamic simulation due to backward 
% Euler integration scheme ('Endpoint').
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
% displacementHooke1DMidpoint
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
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% 
dimension = obj.dimension;
edof = meshObject.edof;
DT = setupObject.timeStepSize;
% material data and voigt notation
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
edRef = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edRef;
uN = edN - edRef;
thetaN1 = dofs.thetaN1;
thetaN = qN(edof(e,:),dimension+1).';
% jacobian
J = qR(edof(e,:),1:dimension)'*dNr';
JN1 = edN1 * dNr';
% initialize energy
elementEnergy.internalEnergy = 0;
elementEnergy.internalEnergyDifference = 0;
%% Gauss loop
for k = 1:numberOfGausspoints
    index = dimension*k-(dimension-1):dimension*k;
    detJ = det(J(:,index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);
    % temperature
    thetaN1e = N(k, :) * thetaN1.';
    thetaNe = N(k, :) * thetaN.';
    % linearized strain
    epsilonN1 = dNx * uN1(:);
    epsilonN = dNx * uN(:);
    % free energy
    PsiMechanicsN1 = 1/2 * epsilonN1 * E * epsilonN1;
    PsiThermoN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiCoupledN1 = -beta*(thetaN1e - thetaR)*trace(epsilonN1);
    PsiN1 = PsiMechanicsN1 + PsiThermoN1 + PsiCoupledN1;
    % entropy
    etaN1 = kappa*log(thetaN1e/thetaR) + beta*epsilonN1;
    etaN = kappa*log(thetaNe/thetaR) + beta*epsilonN;
    % Cauchy stress
    sigmaN1 = E * epsilonN1 - beta*(thetaN1e - thetaR);
    % internal energy
    etaN = kappa*log(thetaNe/thetaR) + beta*epsilonN;
    eN1 = PsiN1 + thetaN1e*etaN1;
    if ~computePostData
        % residual
        RX = RX + dNx' * sigmaN1 * detJ * gaussWeight(k);
        q = -k0*dNx*thetaN1';
        RT = RT + (N(k, :)'*thetaN1e*(etaN1-etaN)/DT-dNx'*q)*detJ*gaussWeight(k);
        % tangent      
        KXX = KXX + dNx'*E*dNx * detJ * gaussWeight(k);
        KXT = KXT - dNx' * beta * N(k,:) * detJ * gaussWeight(k);
        KTX = KTX + N(k,:)' * thetaN1e* 1/DT * beta * dNx * detJ * gaussWeight(k);
        KTT = KTT + (N(k,:)' * ((etaN1-etaN)/DT * N(k,:) + thetaN1e*1/DT*kappa/thetaN1e* N(k,:)) + dNx'*k0*dNx) * detJ * gaussWeight(k);
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + eN1 * detJ * gaussWeight(k);
        elementEnergy.internalEnergyDifference = elementEnergy.internalEnergyDifference + (eN1 - eN) * detJ * gaussWeight(k);
    else
        % stress computation
        stressTensor.Cauchy = sigmaN1;
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
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
