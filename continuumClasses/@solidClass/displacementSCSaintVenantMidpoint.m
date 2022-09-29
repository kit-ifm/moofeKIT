function [rData, kData, elementEnergy, array] = displacementSCSaintVenantMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = stVenantEndPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
% St. Venant Kirchhoff strain-energy function, evaluated at time n+1, i.e. implicid euler method.
%
% 10.01.2012 C.HESCH

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
% mapVoigtObject = obj.mapVoigtObject;
% mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
NAll = shapeFunctionObject.N;
dNrAll = shapeFunctionObject.dNr;
qR = obj.qR;
qN = obj.qN;
numberOfDOFs = size(globalFullEdof, 2);
dimension = obj.dimension;

lambda = materialObject.lambda;
mu = materialObject.mu;
DMat = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
    lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
    lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
    0, 0, 0, mu, 0, 0; ...
    0, 0, 0, 0, mu, 0; ...
    0, 0, 0, 0, 0, mu];
I = eye(dimension);

%% Create residual and tangent
elementEnergy.strainEnergy = 0;

Re = rData{1};
Ke = kData{1, 1};
edN = qN(edof(e, :), 1:dimension)';
edN1 = dofs.edN1;
edN05 = 1 / 2 * (edN + edN1);
J = qR(edof(e, :), 1:dimension)' * dNrAll';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, indx)') \ dNrAll(indx, :);
    FAkt = edN1 * dNx';
    F05 = edN05 * dNx';
    % B-matrix
    B05 = BMatrix(dNx, F05);
    BAkt = BMatrix(dNx, FAkt);
    % Cauchy-Green tensor
    CAkt = FAkt' * FAkt;
    C05 = F05' * F05;
    % Green-Lagrange tensor
    En1 = 0.5 * (CAkt - I);
    En1_v = matrixToVoigt(En1, 'strain');
    En05 = 0.5 * (C05 - I);
    En05_v = matrixToVoigt(En05, 'strain');
    % Strain energy
    elementEnergy.strainEnergy = elementEnergy.strainEnergy + 0.5 * En1_v' * DMat * En1_v * detJ * gaussWeight(k);
    % Stresses
    DW05_v = 1 / 2 * DMat * En05_v;
    DW05 = voigtToMatrix(DW05_v, 'stress');
    % Residual
    Re = Re + 2 * B05' * DW05_v * detJ * gaussWeight(k);
    % Tangent
    D2W05 = 1 / 4 * DMat;
    A1 = 2 * dNx' * DW05 * dNx * detJ * gaussWeight(k);
    MAT = zeros(numberOfDOFs);
    for g = 1:dimension
        MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = 1 / 2 * A1;
    end
    Ke = Ke + 1 / 2 * 4 * B05' * D2W05 * B05 * detJ * gaussWeight(k) + MAT;
end
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end