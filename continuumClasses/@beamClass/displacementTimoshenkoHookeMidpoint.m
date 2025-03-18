function [rData, kData, elementEnergy, array] = displacementTimoshenkoHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
%BEAMTIMOSHENKOENDPOINT Timoshenko beam

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;
shapeFunctionObject = obj.shapeFunctionObject;
mapVoigtObject = obj.mapVoigtObject;

edof = meshObject.edof;
dimension = obj.dimension;
assert(dimension == 1, 'Beam elements are defined only for 1D!');

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

% w and phi at the nodes of the element
wN1 = dofs.edN1;
phiN1 = dofs.phiN1;
edRef = meshObject.nodes(edof(e, :), 1:dimension).';
wN = obj.qN(edof(e, :), 1);
phiN = obj.qN(edof(e, :), 2);
uN = [wN(1);phiN(1);wN(2);phiN(2)];
uN1 = [wN1(1);phiN1(1);wN1(2);phiN1(2)];
uN05 = ( uN + uN1 ) / 2;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

%% Full integration for postprocessing
full_numberOfGausspoints = 2;
[gaussPoints_full, gaussWeight_full] = gaussPointsAndWeights(dimension, full_numberOfGausspoints, 'oneDimensional');
N_k_I_full = computeLagrangeShapeFunction(dimension,2,full_numberOfGausspoints,gaussPoints_full);
         
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

% Jacobian matrices
JAll = computeJacobianForAllGausspoints(edRef, dN_xi_k_I);

% initialize residual & tangent
RW = rData{1};
RPhi = rData{2};
KWW = kData{1, 1};
KWPhi = kData{1, 2};
KPhiW = kData{2, 1};
KPhiPhi = kData{2, 2};

%bending and shear stiffness matrices
stiffness_mat_shear = zeros(4,4); % in GFE: K_s
stiffness_mat_bending = zeros(4,4); % in GFE: K_b

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% Gauss loop
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % B-Matrices
    B_bending = [0, dN_X_I(1), 0, dN_X_I(2)];
    B_shear = [dN_X_I(1), N_k_I(k,1), dN_X_I(2), N_k_I(k,2)];

    if ~computePostData

        %stiffness matrix
        stiffness_mat_shear = stiffness_mat_shear + shearCorFact*G*A*(B_shear')*B_shear * detJ * gaussWeight(k);
        stiffness_mat_bending = stiffness_mat_bending + E*I*(B_bending')*B_bending * detJ * gaussWeight(k); 
        stiffness_mat = stiffness_mat_shear + stiffness_mat_bending;
        
        %residual (linear)
        R = stiffness_mat*uN05;
        
        % strain energy
        kappaN1 = B_bending*uN1; %curvature
        gammaN1 = B_shear*uN1; % shear deformation
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1/2 * (E*I * kappaN1^2 + shearCorFact*G*A*gammaN1^2) * detJ * gaussWeight(k);
        
    end

end

if computePostData

    for k = 1:numberOfGausspoints
        [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
        dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    
        % B-Matrices
        B_bending = [0, dN_X_I(1), 0, dN_X_I(2)];
        B_shear = [dN_X_I(1), N_k_I(k,1), dN_X_I(2), N_k_I(k,2)];
        kappaN1 = B_bending*uN1; %curvature
        gammaN1 = B_shear*uN1; % shear deformation
    end

    % stress resultants
    bending_moment = E*I*kappaN1;
    shear_force = shearCorFact*G*A*gammaN1;

    for k = 1:full_numberOfGausspoints
    
        stressTensor.Cauchy = [bending_moment, shear_force];
        array = postStressComputation(array,N_k_I_full,k,gaussWeight_full,detJ,stressTensor,setupObject,dimension);

    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % requires decomposition into w- and phi-parts
    
    % pass residual
    rData{1} = RW + [R(1);R(3)];
    rData{2} = RPhi + [R(2);R(4)];

    % pass tangent
    kData{1, 1} = KWW + 0.5*[stiffness_mat(1,1), stiffness_mat(1,3); stiffness_mat(3,1), stiffness_mat(3,3)];
    kData{1, 2} = KWPhi + 0.5*[stiffness_mat(1,2), stiffness_mat(1,4); stiffness_mat(3,2), stiffness_mat(3,4)];
    kData{2, 1} = KPhiW + 0.5*[stiffness_mat(2,1), stiffness_mat(2,3); stiffness_mat(4,1), stiffness_mat(4,3)];
    kData{2, 2} = KPhiPhi + 0.5*[stiffness_mat(2,2), stiffness_mat(2,4); stiffness_mat(4,2), stiffness_mat(4,4)];
end

end
