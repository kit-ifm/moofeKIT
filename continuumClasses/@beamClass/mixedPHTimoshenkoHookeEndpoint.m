function [rData, kData, elementData, array] = mixedPHTimoshenkoHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
%BEAMTIMOSHENKOENDPOINT Timoshenko beam

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;
shapeFunctionObject = obj.shapeFunctionObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;
h = setupObject.timeStepSize;

edof = meshObject.edof;
dimension = obj.dimension;
assert(dimension == 1, 'Beam elements are defined only for 1D!');

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

% w and phi at the nodes of the element
wN1 = dofs.edN1;
phiN1 = dofs.phiN1;
edRef = meshObject.nodes(edof(e, :), 1:dimension).';
wN = obj.qN(edof(e, :), 1)';
phiN = obj.qN(edof(e, :), 2)';

% mixed strain quantities
Phi = mixedFEObject.shapeFunctionObject.N_k_I;
gammaN = mixedFEObject.qN(e,1);
kappaN = mixedFEObject.qN(e,2);
gammaN1 = dofs.edAlphaN1(1);
kappaN1 = dofs.edAlphaN1(2);

% material data
E = materialObject.E;
G = materialObject.G;
I = materialObject.I;
A = materialObject.A;
if isfield(materialObject,'shearCorrectionCoefficient')
    shearCorFact = materialObject.shearCorrectionCoefficient;
else
    shearCorFact = 1;
end

% stress resultants
shear_forceN1 = shearCorFact*G*A*gammaN1;
bending_momentN1 = E*I*kappaN1;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% Jacobian matrices
JAll = computeJacobianForAllGausspoints(edRef, dN_xi_k_I);

% initialize elementEnergy
elementData.strainEnergy = 0;

%% Gauss loop
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    L = N_k_I(k, :)' * Phi(k, :);
    K = dN_X_I(:)' * Phi(k, :);
    B = Phi(k, :)' * Phi(k, :);

    if ~computePostData
 
        %residual (linear)
        rData{1} = rData{1} + K'*shear_forceN1 * detJ * gaussWeight(k);
        rData{2} = rData{2} + (K'*bending_momentN1 - L*shear_forceN1) * detJ * gaussWeight(k);
        rData{3} = rData{3} + (B*(gammaN1-gammaN)/h - K*(wN1'-wN')/h + L'*(phiN1'-phiN')/h) * detJ * gaussWeight(k);
        rData{4} = rData{4} + (B*(kappaN1-kappaN)/h - K*(phiN1'-phiN')/h ) * detJ * gaussWeight(k);
        
        %tangent (constant)
        kData{1,3} = kData{1,3} + K'*shearCorFact*G*A * detJ * gaussWeight(k);
        kData{2,3} = kData{2,3} - L *shearCorFact*G*A * detJ * gaussWeight(k);
        kData{2,4} = kData{2,4} + K'*E*I * detJ * gaussWeight(k);
        kData{3,3} = kData{3,3} + 1/h*B * detJ * gaussWeight(k);
        kData{4,4} = kData{4,4} + 1/h*B * detJ * gaussWeight(k);
        kData{3,1} = kData{3,1} - 1/h*K * detJ * gaussWeight(k);
        kData{3,2} = kData{3,2} + 1/h*L'* detJ * gaussWeight(k);
        kData{4,2} = kData{4,2} - 1/h*K * detJ * gaussWeight(k);
        
        % % strain energy
        elementData.strainEnergy = elementData.strainEnergy + 1/2 * (E*I * kappaN1^2 + shearCorFact*G*A*gammaN1^2) * detJ * gaussWeight(k);
        
    else
        
        stressTensor.Cauchy = [bending_momentN1, shear_forceN1];
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);

    end

end

end