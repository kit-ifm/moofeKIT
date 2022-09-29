function [rData, kData, elementEnergy, array] = displacementHookeSplitMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
%% Creates the residual and the tangent of the given obj.
% 08.02.2015 Marlon Franke: first esra version
% 19.08.2021 Marlon Franke: moofeKIT version

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
% mixedFEObject = obj.mixedFEObject;

meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
dimension = obj.dimension;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% nodal dofs
qR = obj.qR;
qN = obj.qN;
qN1 = dofs.edN1;

% Visco Information
DT = setupObject.timeStepSize;
eModul0 = materialObject.eModul0;
eModul1 = materialObject.eModul1;
eta1 = materialObject.eta1;
tau = eta1/eModul1;
%innere Variable stressVisco hier auch als epsilonVisco abgespeichert
epsilonViscoN = obj.epsilonViscoN;
epsilonViscoN1 = obj.epsilonViscoN1;

eModul = eModul1 + eModul0;
gamma = eModul1/eModul;

% material data and voigt notation
lambda = materialObject.lambda;
mu = materialObject.mu;
selectMapVoigt(mapVoigtObject,dimension,'symmetric');

if (dimension == 1)
    if isfield(materialObject,'E')
        C = materialObject.E;
    else
        C = eModul - (DT*gamma)/(2*tau+DT) * eModul;
    end
else
    C = [lambda+2*mu lambda lambda 0 0 0;...
        lambda lambda+2*mu lambda 0 0 0;...
        lambda lambda lambda+2*mu 0 0 0;...
        0 0 0 mu 0 0;...
        0 0 0 0 mu 0;...
        0 0 0 0 0 mu];
end


edN = qN(edof(e, :), 1:dimension)';
edN1 = qN1;
edR = qR(edof(e, :), 1:dimension)';
uN = edN - edR;
uN1 = edN1 - edR;
uN05 = 1 / 2 * (uN1 + uN);
stressViscoN = epsilonViscoN(e,:);
stressViscoN1 = zeros(1,numberOfGausspoints);
stressViscoN05 = zeros(1,numberOfGausspoints);

% initialize energy
elementEnergy.strainEnergy = 0;
elementEnergy.dissipatedWork = 0;

%% Initialization
RX = rData{1, 1};
KXX = kData{1, 1};

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);    
    B = BMatrix(dN_X_I,'mapVoigtObject',mapVoigtObject);
    if ~computePostData
        epsilonN05 = B*uN05(:);
        epsilonN1 = B*uN1(:);
% Newton-Schleife für innere Variable stressVisco
newtonStep = 0;
stressViscoN1(k) = stressViscoN(k);
    while newtonStep <= setupObject.newton.maximumSteps
        newtonStep = newtonStep + 1; 
        R=DT/tau*gamma*eModul*epsilonN05 - (1+DT/(2*tau))*stressViscoN1(k)+(1-DT/(2*tau))*stressViscoN(k);
        K=-1-DT/(2*tau);
        normR = norm(R);
%         fprintf('|      iteration %2d: norm(R) = %9.5e    |\n',newtonStep,normR)
        if (newtonStep > setupObject.newton.maximumSteps || normR > setupObject.newton.maximumResiduumNorm)
            error('Newton iteration failed!');
        elseif normR < setupObject.newton.tolerance
            % converged state
            break;
         else
            % direct solver
            stressViscoN1(k) = stressViscoN(k) - R/K;
        end
    end
%         stressViscoN1(k) = (DT/tau*gamma*eModul*epsilonN05+(1-DT/(2*tau))*stressViscoN(k))/(1+DT/(2*tau));
        epsilonViscoN1(e,k) = stressViscoN1(k);
        stressViscoN05(k) = 1/2 *(stressViscoN(k)+stressViscoN1(k));
        % dissipative Arbeit
        elementEnergy.dissipatedWork = elementEnergy.dissipatedWork + DT*(epsilonN05-stressViscoN05(k)/eModul1)*eModul1*eModul1*detJ*gaussWeight(k)*(epsilonN05-stressViscoN05(k)/eModul1)/eta1;
        % Spannung
            %Anfangsspannungen: sigma0
            sigma0N05 = eModul*epsilonN05;
        sigmaN05 = sigma0N05 - 1/2 *(stressViscoN(k)+stressViscoN1(k));
        % Residual
        RX = RX + B'.*sigmaN05*detJ*gaussWeight(k);
        % Tangent
        KXX = KXX + 0.5*B'*C*B*detJ*gaussWeight(k);
        %strain energy 
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1/2*(epsilonN1*eModul0*epsilonN1 + (epsilonN1-stressViscoN1(k)/eModul1)*eModul1*(epsilonN1-stressViscoN1(k)/eModul1))*detJ*gaussWeight(k);
    else
        % stress at gausspoint
        epsilon = B*uN05(:);
        sigma_V = C*epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJStruct,stressTensor,setupObject,dimension);
    end
end
obj.epsilonViscoN1 = epsilonViscoN1;
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end