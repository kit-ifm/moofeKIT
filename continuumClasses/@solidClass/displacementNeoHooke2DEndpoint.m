function [rData, kData, elementEnergy, array] = displacementNeoHooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
% material data and voigt notation
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
I = eye(dimension);
E = materialObject.E;
nu = materialObject.nu;
mu = E / (2 * (1 + nu));
beta = 2 * nu / (1 - 2 * nu);
% initialization residual, tangent & elementEnergy
RX = rData{1};
KXX = kData{1, 1};
elementEnergy.strainEnergy = 0;
%loop over all gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    %deformation gradient
    FN1 = edN1 * dN_X_I';
    detFN1 = det(FN1);

    %nodal operatormatix
    B = BMatrix(dN_X_I, FN1, 'mapVoigtObject', mapVoigtObject);

    %right cauchy-green
    CN1 = FN1' * FN1;
    CN1inv = CN1 \ I;
    CN1invSym = [CN1inv .* CN1inv, CN1inv(:, 1) .* CN1inv(:, 2); CN1inv(1, :) .* CN1inv(2, :), (CN1inv(1, 1) * CN1inv(2, 2) + CN1inv(2, 1) * CN1inv(1, 2)) / 2];
    CN1invVoigt = [CN1inv(1, 1), CN1inv(2, 2), CN1inv(2, 1)]';

    %derivatives of strain energy
    psi = mu / 2 * (trace(CN1) - dimension) + mu / beta * (detFN1^-beta - 1);
    psi_C = mu / 2 * (I - detFN1^(-beta) * CN1inv);
    psi_CC = mu / 2 * (beta / 2 * detFN1^(-beta) * (CN1invVoigt * CN1invVoigt') + detFN1^(-beta) * CN1invSym);
    %2.PK
    SN1 = 2 * psi_C;
    SN1Voigt = [SN1(1, 1); SN1(2, 2); SN1(1, 2)];
    if ~computePostData
        %tangent and rediduum
        KXX = KXX + (B' * 4 * psi_CC * B + kron(dN_X_I'*SN1*dN_X_I, I)) * detJ * gaussWeight(k);
        %         Re = Re + (B'*SN1Voigt-Ni'*b) * detJ*gaussWeight(j);
        RX = RX + (B' * SN1Voigt) * detJ * gaussWeight(k);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + psi * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end
