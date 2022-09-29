function [rData, kData, elementEnergy, array] = displacementSCSaintVenantEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCSAINTVENANTENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic  Saint-Venant
% Kirchhoff ('SaintVenant') model (attention: nonlinear geometric but
% linear material/stress-strain relation).
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementSCSaintVenantEndpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
%
% REFERENCE
% Nichtlineare Finite-Element-Methoden, Peter Wriggers, Springer-Verlag,
% 2001, https://link.springer.com/book/10.1007/978-3-642-56865-7
%
% SEE ALSO
% displacementSCSaintVenantMidpoint,
% displacementSCSaintVenantAutomaticDiffEndpoint,
% displacementSCSaintVenantAutoDiffEndpoint
%
% CREATOR(S)
% Marlon Franke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
% mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfDOFs = size(globalFullEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
% material data and voigt notation
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
lambda = materialObject.lambda;
mu = materialObject.mu;
DMat = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
    lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
    lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
    0, 0, 0, mu, 0, 0; ...
    0, 0, 0, 0, mu, 0; ...
    0, 0, 0, 0, 0, mu];
I = eye(3);

%% Create residual and tangent
elementEnergy.strainEnergy = 0;
Re = rData{1};
Ke = kData{1, 1};

% post processing stresses
edN1 = dofs.edN1;
J = qR(edof(e, :), 1:dimension).' * dNr';
JN1 = edN1 * dNr';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    index = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, index).');
    detJN1 = det(JN1(:, index).');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index).') \ dNr(index, :);
    FN1 = edN1 * dNx';
    % B-matrix
    BN1 = BMatrix(dNx, FN1);
    % Cauchy-Green tensor
    CN1 = FN1.' * FN1;
    % Green-Lagrange tensor
    EN1 = 0.5 * (CN1 - I);
    EN1_v = matrixToVoigt(EN1, 'strain');
    % Stresses
    DW1_v = 1 / 2 * DMat * EN1_v;
    DW1 = voigtToMatrix(DW1_v, 'stress');
    if ~computePostData
        % Residual
        Re = Re + 2 * BN1' * DW1_v * detJ * gaussWeight(k);
        % Tangent
        D2W1 = 1 / 4 * DMat;
        A1 = 2 * dNx' * DW1 * dNx * detJ * gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        Ke = Ke + 4 * BN1' * D2W1 * BN1 * detJ * gaussWeight(k) + MAT;
        % Strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 0.5 * EN1_v' * DMat * EN1_v * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        SN1_v = DMat * EN1_v;
        SN1 = voigtToMatrix(SN1_v, 'stress');
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end