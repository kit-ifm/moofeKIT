function [rData, kData, elementEnergy, array] = displacementHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% homogenous linear-elastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
% 08.02.2015 Marlon Franke: first esra version
% 19.08.2021 Marlon Franke: moofeKIT version
%creates the residual and the tangent of the given obj.
%

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
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;

% nodal dofs
qR = obj.qR;
qN1 = dofs.edN1;

% material data and voigt notation
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

edN1 = qN1;
edRef = qR(edof(e, :), 1:dimension)';
uN1 = edN1 - edRef;
J = qR(edof(e, :), 1:dimension)' * dNr';
JN1 = qN1 * dNr';

% initialize energy
elementEnergy.strainEnergy = 0;

%% Initialization
RX = rData{1, 1};
KXX = kData{1, 1};

%% Gauss loop
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    detJN1 = det(JN1(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, indx)') \ dNr(indx, :);
    B = BMatrix(dNx, 'mapVoigtObject', mapVoigtObject);
    if ~computePostData
        % Residual
        RX = RX + (B' * C * B) * uN1(:) * detJ * gaussWeight(k);
        % Tangent
        KXX = KXX + B' * C * B * detJ * gaussWeight(k);
        %strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1(:))' * C * (B * uN1(:)) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        epsilon = B * uN1(:);
        sigma_V = C * epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end