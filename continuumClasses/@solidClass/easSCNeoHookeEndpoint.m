function [rData, kData, elementEnergy, array] = easSCNeoHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = neoHookEndPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
% Neo Hook strain-energy function, enhanced modes, evaluated at time n+1, i.e. implicid euler method.
%
% 10.01.2012 C.HESCH
% 25.8.2018 Resolved issue with tangent and added 2D routine R.PFEFFERKORN

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
mapVoigtObject = obj.mapVoigtObject;
meshObject = obj.meshObject;

% aquire general data
N = shapeFunctionObject.N;
dNrAll = shapeFunctionObject.dNr;
dNr0 = shapeFunctionObject.dNr0;
dMr = mixedFEObject.shapeFunctionObject.dNr;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;

dimension = obj.dimension;

I = eye(3);

% aquire material data
lambda = materialObject.lambda;
mu = materialObject.mu;

% aquire the nodal values of the variables for the current element
edRef = obj.qR(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
edAlphaN1 = dofs.edAlphaN1;

numberOfDisplacementDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalDofs = size(mixedFEObject.qN1, 2);

if dimension == 2
    alphaN1Voigt = [edAlphaN1(1, 1:2:end); edAlphaN1(1, 2:2:end)];
else
    alphaN1Voigt = [edAlphaN1(1, 1:3:end); edAlphaN1(1, 2:3:end); edAlphaN1(1, 3:3:end)];
end

% compute Jacobian matrices
J = edRef * dNrAll';
JN1 = edN1 * dNrAll';
J0 = edRef * dNr0';

% initialize residual
RD = rData{1};
RA = rData{2};

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    detJN1 = det(JN1(:, indx)');
    detJ0 = det(J0);
    dNx = (J(:, indx)') \ dNrAll(indx, :);
    dMx = (detJ0 / detJ) * (J0' \ dMr(indx, :));
    % deformation gradient
    FAkt = eye(3);
    FAkt(1:dimension, 1:dimension) = edN1 * dNx' + alphaN1Voigt * dMx';
    % B-matrix (current configuration)
    BAkt = BMatrix(dNx, FAkt);
    GAkt = BMatrix(dMx, FAkt);
    % Cauchy-Green tensor
    CAkt = FAkt' * FAkt;
    % Inverse strain tensor
    CInvAkt = CAkt \ I;
    CInvAktSym = [CInvAkt .* CInvAkt, CInvAkt .* CInvAkt(:, [2, 3, 1]); CInvAkt .* CInvAkt([2, 3, 1], :), (CInvAkt .* CInvAkt([2, 3, 1], [2, 3, 1]) + CInvAkt(:, [2, 3, 1]) .* CInvAkt([2, 3, 1], :)) / 2];
    CInvAkt_v = matrixToVoigt(CInvAkt, 'stress');
    % Invariant
    IAkt = (det(CAkt))^0.5;
    % Strain energy function
    elementEnergy.strainEnergy = elementEnergy.strainEnergy + (mu / 2 * (trace(CAkt) - 3) + lambda / 2 * (log(IAkt))^2 - mu * log(IAkt)) * detJ * gaussWeight(k);
    % First derivative of strain energy function
    DWAkt = mu / 2 * (I - CInvAkt) + lambda / 2 * log(IAkt) * CInvAkt;
    DWAkt = DWAkt(1:dimension, 1:dimension);
    DWAkt_v = matrixToVoigt(DWAkt, 'stress');
    % Second derivative of strain energy function
    D2W1 = (mu / 2 - lambda / 2 * log(IAkt)) * CInvAktSym + lambda / 4 * (CInvAkt_v * CInvAkt_v');
    if dimension == 2
        D2W1 = D2W1([1, 2, 4], [1, 2, 4]);
    end
    if ~computePostData
        % RESIDUAL
        RD = RD + 2 * BAkt' * DWAkt_v * detJ * gaussWeight(k);
        RA = RA + 2 * GAkt' * DWAkt_v * detJ * gaussWeight(k);

        % TANGENT
        A1 = 2 * dNx' * DWAkt * dNx * detJ * gaussWeight(k);
        A2 = 2 * dNx' * DWAkt * dMx * detJ * gaussWeight(k);
        A3 = 2 * dMx' * DWAkt * dMx * detJ * gaussWeight(k);
        MAT1 = zeros(numberOfDisplacementDofs);
        MAT2 = zeros(numberOfDisplacementDofs, numberOfInternalDofs);
        MAT3 = zeros(numberOfInternalDofs);
        for g = 1:dimension
            MAT1(g:dimension:numberOfDisplacementDofs, g:dimension:numberOfDisplacementDofs) = A1;
        end
        for g = 1:dimension
            MAT2(g:dimension:numberOfDisplacementDofs, g:dimension:numberOfInternalDofs) = A2;
        end
        for g = 1:dimension
            MAT3(g:dimension:numberOfInternalDofs, g:dimension:numberOfInternalDofs) = A3;
        end
        KDD = KDD + 4 * BAkt' * D2W1 * BAkt * detJ * gaussWeight(k) + MAT1;
        KDA = KDA + 4 * BAkt' * D2W1 * GAkt * detJ * gaussWeight(k) + MAT2;
        KAA = KAA + 4 * GAkt' * D2W1 * GAkt * detJ * gaussWeight(k) + MAT3;
    else
        % STRESS COMPUTATION
        PN1 = 2 * FAkt * DWAkt;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / (det(FAkt)) * PN1 * FAkt';
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RD;
    rData{2} = RA;

    % pass tangent
    kData{1, 1} = KDD;
    kData{1, 2} = KDA;
    kData{2, 1} = KDA';
    kData{2, 2} = KAA;
end
end