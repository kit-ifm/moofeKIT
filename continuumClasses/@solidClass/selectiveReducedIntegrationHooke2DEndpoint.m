function [rData, kData, elementEnergy, array] = selectiveReducedIntegrationHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
numberOfDofs = size(globalFullEdof, 2);
dimension = obj.dimension;
% nodal dofs
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
% material data
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
nu = materialObject.nu;
E = materialObject.E;
mu = E / (2 * (1 + nu));
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
Iv = [1, 1, 0]';
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'evz'
        Cmu = mu * diag([2, 2, 1]);
        Clam = lambda * (Iv * Iv');
    otherwise
        error('not implemented')
end
%nodal points and displacements
u = zeros(numberOfDofs, 1);
u(:) = edN1 - edR;

I = eye(dimension);
% initialization residual, tangent & elementEnergy
RX = rData{1};
KXX = kData{1, 1};
elementEnergy.strainEnergy = 0;
%% loop of integration types (reduced and full)
for integrationType = 1:2
    %shape functions
    if integrationType == 1
        numberOfGausspoints = size(mixedFEObject.shapeFunctionObject.gaussWeight, 2);
        gaussWeight = mixedFEObject.shapeFunctionObject.gaussWeight;
        N_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
        dN_xi_k_I = mixedFEObject.shapeFunctionObject.dN_xi_k_I;
    else
        numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
        gaussWeight = shapeFunctionObject.gaussWeight;
        N_k_I = shapeFunctionObject.N_k_I;
        dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
    end
    if (computePostData && integrationType == 2) || ~computePostData
        %loop over all gauss points
        for k = 1:numberOfGausspoints
            %approximation of geometry and deriviates according to x
            [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
            %nodal operatormatix
            B = BMatrix(dN_X_I, 'mapVoigtObject', mapVoigtObject);
            if ~computePostData
                %tangent and residuum
                if integrationType == 1
                    KXX = KXX + B' * Clam * B * detJ * gaussWeight(k);
                    RX = RX + (B' * Clam * B * u) * detJ * gaussWeight(k);
                else
                    KXX = KXX + B' * Cmu * B * detJ * gaussWeight(k);
                    %                 Re = Re + (B'*Cmu*B*u-Ni'*b) * detJ*gaussWeight(j);
                    RX = RX + (B' * Cmu * B * u) * detJ * gaussWeight(k);
                    elementEnergy.strainEnergy = elementEnergy.strainEnergy + (B * u)' * Clam * B * u * detJ * gaussWeight(k) + (B * u)' * Cmu * B * u * detJ * gaussWeight(k);
                end
            else
                gradU = (edN1 - edR) * dN_X_I';
                epsilon = zeros(3, 3);
                epsilon(1:2, 1:2) = 1 / 2 * (gradU + gradU');
                switch lower(strtok(obj.materialObject.name, 'Hooke'))
                    case 'esz'
                        epsilon(3, 3) = -nu / (1 - nu) * (trace(epsilon)); %ESZ
                    case 'evz'
                        epsilon(3, 3) = 0; %EVZ
                    otherwise
                        error('not implemented')
                end
                %stresses
                sigmaN1 = lambda * trace(epsilon) * eye(3) + 2 * mu * epsilon;
                stressTensor.Cauchy = sigmaN1;
                array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
            end
        end
        if ~computePostData
            rData{1} = RX;
            kData{1, 1} = KXX;
        end
    end
end

end