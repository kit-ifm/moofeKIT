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
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% nodal dofs
edR = obj.qR(edof(e, :), :).';
edN = obj.qN(edof(e, :), :).';
edN1 = dofs.edN1;
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

JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

%% Create residual and tangent
elementEnergy.strainEnergy = 0;
Re = rData{1};
Ke = kData{1, 1};

% post processing stresses
if computePostData
    JNAll = computeJacobianForAllGausspoints(edN, dN_xi_k_I);
    JN1All = computeJacobianForAllGausspoints(edN1, dN_xi_k_I);
end
% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    FN1 = edN1 * dN_X_I';
    % B-matrix
    BN1 = BMatrix(dN_X_I, FN1);
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
        A1 = 2 * dN_X_I' * DW1 * dN_X_I * detJ * gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        Ke = Ke + 4 * BN1' * D2W1 * BN1 * detJ * gaussWeight(k) + MAT;
        % Strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 0.5 * EN1_v' * DMat * EN1_v * detJ * gaussWeight(k);
    else
        [~, detJN] = extractJacobianForGausspoint(JNAll, k, setupObject, dimension);
        [~, detJN1] = extractJacobianForGausspoint(JN1All, k, setupObject, dimension);
        detJStruct = struct('R', detJ, 'N', detJN, 'N1', detJN1);
        % stress at gausspoint
        SN1_v = DMat * EN1_v;
        SN1 = voigtToMatrix(SN1_v, 'stress');
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end