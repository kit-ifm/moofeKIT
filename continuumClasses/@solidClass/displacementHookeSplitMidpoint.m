function [rData, kData, elementEnergy, array] = displacementHookeSplitMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
% 08.02.2015 Marlon Franke: first esra version
% 19.08.2021 Marlon Franke: moofeKIT version

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;
plotObject = setupObject.plotObject;

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
qN1 = obj.qN1;
% material data and voigt notation
lambda = materialObject.lambda;
mu = materialObject.mu;
K = lambda + 2/3 * mu;
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
if (dimension == 1)
    if isfield(materialObject, 'E')
        C = materialObject.E;
    else
        C = mu * (3 * lambda + 2 * mu) / (lambda + mu);
    end
else
%     C = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
%         lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
%         lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
%         0, 0, 0, mu, 0, 0; ...
%         0, 0, 0, 0, mu, 0; ...
%         0, 0, 0, 0, 0, mu];
    % neues C definieren: C = C_dev + C_vol
    C_dev = [4/3 * mu, -2/3 * mu, -2/3 * mu, 0, 0, 0; ...
             -2/3 * mu, 4/3 * mu, -2/3 * mu, 0, 0, 0; ...
             -2/3 * mu, -2/3 * mu, 4/3 * mu, 0, 0, 0; ...
             0, 0, 0, mu, 0, 0; ...
             0, 0, 0, 0, mu, 0; ...
             0, 0, 0, 0, 0, mu];
    C_vol = [K, K, K, 0, 0, 0; ...
             K, K, K, 0, 0, 0; ...
             K, K, K, 0, 0, 0; ...
             0, 0, 0, 0, 0, 0; ...
             0, 0, 0, 0, 0, 0; ...
             0, 0, 0, 0, 0, 0];
    C = C_dev + C_vol;
end

%% Create residual and tangent
elementEnergy.strainEnergy = 0;

Re = rData{1};
Ke = kData{1, 1};
edN1 = dofs.edN1;
edN = qN(edof(e, :), 1:dimension)';
edR = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edR;
uN = edN - edR;
uN05 = 1 / 2 * (uN1 + uN);
% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
    epsilon = B*uN05(:);
    IV = [1, 1, 1, 0, 0, 0]; %Einheitsmatrix in Voigt-Notation
    epsilonVol = 1/3*(epsilon(1,1) + epsilon(2,1) + epsilon(3,1))*IV';
    epsilonDev = epsilon - epsilonVol;
    if ~computePostData
        sigmaDev = C_dev*epsilonDev;
        sigmaVol = C_vol*epsilonVol;
        sigma_V = sigmaVol + sigmaDev;
        Re = Re + B' * sigma_V * detJ * gaussWeight(k);
        Ke = Ke + 1 / 2 * B' * C * B * detJ * gaussWeight(k);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1(:))' * C * (B * uN1(:)) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        epsilon = B*uN05(:);
        sigma_V = C*epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);
    end
end
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end