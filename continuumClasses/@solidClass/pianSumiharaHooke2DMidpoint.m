function [rData, kData, elementEnergy, array] = pianSumiharaHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% Pian-Sumihara element ansatz

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% general data
dimension = obj.dimension;
% data for displacement dofs
NAll = shapeFunctionObject.N;
dNrAll = shapeFunctionObject.dNr;
dNr0 = shapeFunctionObject.dNr0;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
numberOfDisplacementDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
qR = obj.qR;
qN = obj.qN;
edof = meshObject.edof;
% mixed FE data
Pr = mixedFEObject.shapeFunctionObject.N;
numberOfInternalDofs = size(mixedFEObject.qN1, 2);
betaN = mixedFEObject.qN;
% material data
nu = materialObject.nu;
E = materialObject.E;
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2]; %EVZ
    otherwise
        error('not implemented')
end
Cinv = C \ eye(size(C, 1));

elementEnergy.strainEnergy = 0;

Ge = zeros(numberOfInternalDofs, numberOfDisplacementDofs);
He = zeros(numberOfInternalDofs, numberOfInternalDofs);
edN = qN(edof(e, :), 1:dimension)';
edN1 = dofs.edN1;
edRef = qR(edof(e, :), 1:dimension)';
uN = edN - edRef;
uN1 = edN1 - edRef;
uN05 = 1 / 2 * (uN + uN1);
betaNe = betaN(e, :)';
betaN1e = dofs.edAlphaN1';
betaN05e = 1 / 2 * (betaN1e+betaNe);

%jacobian element centroid
J0 = edRef * dNr0';
F0 = F0Matrix(dimension, J0);

J = qR(edof(e, :), 1:dimension)' * dNrAll';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, indx)') \ dNrAll(indx, :);

    %nodal operatormatix
    B = BMatrix(dNx, 'mapVoigtObject', mapVoigtObject);

    %stress approximation
    indx = 3 * (k - 1) + 1:3 * k;
    P = F0 * Pr(indx, :);

    %tangent and rediduum
    Ge = Ge + P' * B * detJ * gaussWeight(k);
    He = He - P' * Cinv * P * detJ * gaussWeight(k);

    sigmaVoigt = P * betaN1e;
    epsilonVoigt = Cinv * sigmaVoigt;
    elementEnergy.strainEnergy = elementEnergy.strainEnergy + (sigmaVoigt' * epsilonVoigt - 1 / 2 * sigmaVoigt' * Cinv * sigmaVoigt) * detJ * gaussWeight(k);
end
if ~computePostData
    rData{1} = Ge' * betaN05e;
    rData{2} = Ge * uN05(:) + He * betaN05e;
    kData{1, 1} = zeros(numberOfDisplacementDofs);
    kData{1, 2} = 0.5 * Ge';
    kData{2, 1} = 0.5 * Ge;
    kData{2, 2} = 0.5 * He;
end
end
