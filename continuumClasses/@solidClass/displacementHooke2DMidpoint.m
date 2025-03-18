function [rData, kData, elementEnergy, array] = displacementHooke2DMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
% nodal dofs
edR = obj.qR(edof(e, :), 1:dimension).';
edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
% material data
nu = materialObject.nu;
E = materialObject.E;
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'esz'
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; %ESZ
    case 'evz'
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2]; %EVZ
    otherwise
        error('not implemented')
end

% Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

%% Create residual and tangent
elementEnergy.strainEnergy = 0;

Re = rData{1};
Ke = kData{1, 1};
uN = edN - edR;
uN1 = edN1 - edR;
uN05 = 1 / 2 * (uN + uN1);
% Run through all Gauss points
for k = 1:numberOfGausspoints
    % compute the Jacobian determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    B = BMatrix(dN_X_I);
    if ~computePostData
        Re = Re + (B' * C * B) * uN05(:) * detJ * gaussWeight(k);
        Ke = Ke + 1 / 2 * B' * C * B * detJ * gaussWeight(k);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1(:))' * C * (B * uN1(:)) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        [detJ, detJStruct, dN_X_k_I, ~] = computeAllJacobian(edR, edN, edN1, dN_xi_k_I, k, setupObject);

        sigmaN1_v = C * B * uN1(:);
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end