function [rData, kData, elementEnergy, array] = displacementANSTimoshenkoHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
% DISPLACEMENTANSTIMOSHENKOENDPOINT Element routine of class beamClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% In this formulation, the transverse shear locking is eliminated via an
% Assumed Natural Strain (ANS) approach.
%
% CALL
% displacementPetrovGalerkinANSTimoshenkoEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type beamClass,
%      e.g. beamObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data.
% kData: cell-array of size [2, 2] for tangent data.
% dofs: degrees of freedom (dofs) optionally manipulated data (numerical
%       tangent)
% array: structure for storage fe-data, for more information see
%        storageFEObject.initializeArrayStress
% stressTensor: structure for storage stress tensors (postprocessing), for
%               more information see storageFEObject.initializeArrayStress
% flagNumericalTangent: flag that indicates whether the function call
%                       happens during the computation of the numerical
%                       tangent or not.
%
% REFERENCE
% -
%
% SEE ALSO
% -
%
% CREATOR(S)
% Felix Zaehringer

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;
shapeFunctionObject = obj.shapeFunctionObject;

edof = meshObject.edof;
dimension = obj.dimension;
assert(dimension == 1, 'Beam elements are defined only for 1D!');

% w and phi at the nodes of the element
wN1 = dofs.edN1;
phiN1 = dofs.phiN1;
edRef = meshObject.nodes(edof(e, :), :).';
uN1 = [wN1(1);phiN1(1);wN1(2);phiN1(2)];

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
         
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

he = norm(edRef(:, 1)-edRef(:, 2));

% initialize residual & tangent
RW = rData{1};
RPhi = rData{2};
KWW = kData{1, 1};
KWPhi = kData{1, 2};
KPhiW = kData{2, 1};
KPhiPhi = kData{2, 2};

%bending and shear stiffness matrices
stiffness_mat_shear = zeros(4, 4); % in GFE: K_s
stiffness_mat_bending = zeros(4, 4); % in GFE: K_b

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% Gauss loop
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % B-Matrices
    B_bending = [0, dN_X_I(1), 0, dN_X_I(2)];
    B_shear = [-1/he, 1/2, 1/he, 1/2];

    if ~computePostData

        %stiffness matrix
        stiffness_mat_shear = stiffness_mat_shear + shearCorFact*G*A*(B_shear')*B_shear * detJ * gaussWeight(k);
        stiffness_mat_bending = stiffness_mat_bending + E*I*(B_bending')*B_bending * detJ * gaussWeight(k); 
        stiffness_mat = stiffness_mat_shear + stiffness_mat_bending;
        
        %residual (linear)
        R = stiffness_mat*uN1;
        
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

    for k = 1:numberOfGausspoints
    
        stressTensor.Cauchy = [bending_moment, shear_force];
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);

    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % requires decomposition into w- and phi-parts
    
    % pass residual
    rData{1} = RW + [R(1);R(3)];
    rData{2} = RPhi + [R(2);R(4)];

    % pass tangent
    kData{1, 1} = KWW + stiffness_mat(1:2:end, 1:2:end);
    kData{1, 2} = KWPhi + stiffness_mat(1:2:end, 2:2:end);
    kData{2, 1} = KPhiW + stiffness_mat(2:2:end, 1:2:end);
    kData{2, 2} = KPhiPhi + stiffness_mat(2:2:end, 2:2:end);
end

end

