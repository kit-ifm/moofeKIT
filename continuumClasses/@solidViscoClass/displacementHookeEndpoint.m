function [rData, kData, elementEnergy, array] = displacementHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
epsilonViscoN = obj.epsilonViscoN;
epsilonViscoN1 = obj.epsilonViscoN1;
eVN1 = zeros(1,numberOfGausspoints);
% material data and voigt notation
lambda = materialObject.lambda;
mu = materialObject.mu;
selectMapVoigt(mapVoigtObject,dimension,'symmetric');

if (dimension == 1)
    if isfield(materialObject,'E')
        C = materialObject.E;
    else
        C = eModul0 + eModul1/(1+DT/(2*tau));
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
eVN = epsilonViscoN(e,:);

% initialize energy
elementEnergy.strainEnergy = 0;

%% Initialization
RX = rData{1, 1};
KXX = kData{1, 1};


% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    B = BMatrix(dN_X_I,'mapVoigtObject',mapVoigtObject);
    if ~computePostData
        epsilon = B*uN1(:);
        % innere Variable epsilonVisco
        eVN1(k) = (eVN(k)*(1-DT*eModul1/(2*eta1))+((DT*eModul1)/eta1)*epsilon)/(1+DT*eModul1/(2*eta1));
        epsilonViscoN1(e,k) = eVN1(k);
        % Spannung
        sigma = (eModul0+eModul1/(1+DT/(2*tau)))*epsilon-eModul1/(1+DT/(2*tau))*eVN(k);
        % Residual
        RX = RX + B'.*sigma*detJ*gaussWeight(k);
        % Tangent
        KXX = KXX + B'*C*B*detJ*gaussWeight(k);
        %strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1/2*(B*uN1(:))'*C*(B*uN1(:))*detJ*gaussWeight(k);
    else
        % stress at gausspoint
        epsilon = B*uN1(:);
        sigma_V = C*epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);
    end
end
obj.epsilonViscoN1 = epsilonViscoN1;
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end