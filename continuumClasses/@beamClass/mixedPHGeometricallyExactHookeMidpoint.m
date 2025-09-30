function [rData, kData, elementData, array] = mixedPHGeometricallyExactHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
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

% r and phi (displacement quantities)
edN1  = dofs.edN1;
rN1   = edN1(:);
phiN1 = dofs.phiN1';
edN   = obj.qN(edof(e, :), 1:2)';
rN    = edN(:);

%%%%
dimAbs = 2;
edR = obj.qR(edof(e, :), 1:dimAbs).';

%%%%%

phiN   = obj.qN(edof(e, :), 3);
rN05   = 1/2*(rN + rN1);
phiN05 = 1/2*(phiN + phiN1);

% mixed strain quantities
Psi = mixedFEObject.shapeFunctionObject.N_k_I;
GammaN = mixedFEObject.qN(e,1:2)';
kappaN = mixedFEObject.qN(e,3);
GammaN1 = dofs.edAlphaN1(1:2)';
kappaN1 = dofs.edAlphaN1(3);
kappaN05 = 1/2*(kappaN + kappaN1);
GammaN05 = 1/2*(GammaN + GammaN1);

% material data
[EA, EI, kGA] = get_material_parameters(materialObject);

% matrices
skew_matrix = [0 , 1;
    -1, 0];
material_matrix = [EA, 0;
    0, kGA];

% stress resultants
resultant_forcesN05 = material_matrix*GammaN05;
resultant_momentN05 = EI*kappaN05;
resultant_forcesN1  = material_matrix*GammaN1;
resultant_momentN1  = EI*kappaN1;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% initialize elementEnergy
elementData.strainEnergy = 0;

% coefficient matrices
G1  = zeros(2,1);
GR2 = zeros(4,2);
C1  = zeros(1,1);
C2  = zeros(2,2);
V   = zeros(2,2);
%M2  = zeros(4,4);
%M1  = zeros(2,2);

%% Gauss loop
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I    = computedN_X_I(dN_xi_k_I, J, k);
    
    % shape function matrices
    Psi1   = Psi(k,:);
    Psi2   = [Psi(k,:), 0; 0, Psi(k,:)];
    Phi1   = N_k_I(k, :);
    dsPhi1 = dN_X_I;
    dsPhi2 = [dN_X_I(1) 0         dN_X_I(2) 0        ;
        0         dN_X_I(1) 0         dN_X_I(2)];
    %Phi2 = [Phi1(1) 0         Phi1(2) 0        ;
    %          0         Phi1(1) 0         Phi1(2)];
    
    
    if ~computePostData
        % rotation matrix
        RN05 = get_planar_rotation_matrix(Phi1*phiN05);
        RN1 = get_planar_rotation_matrix(Phi1*phiN1);
        
        % coefficient matrices
        G1  = G1  + dsPhi1' * Psi1 * detJ * gaussWeight(k);
        C1  = C1  + Psi1'   * Psi1 * detJ * gaussWeight(k);
        C2  = C2  + Psi2'   * Psi2 * detJ * gaussWeight(k);
        GR2 = GR2 + dsPhi2' * RN05 * Psi2 * detJ * gaussWeight(k) ;
        V   = V   + Phi1'   * (dsPhi2*rN05)' * skew_matrix * RN05 * Psi2 * detJ * gaussWeight(k);
        
        % potential energy due to strain
        elementData.strainEnergy = elementData.strainEnergy + 1/2 * ((GammaN1)'*material_matrix*(GammaN1) + EI * kappaN1^2) * detJ * gaussWeight(k);
        elementData.drift = norm(GammaN1 - RN1'*(dsPhi2*rN1) +[1;0]);
    else
        
        stressTensor.Cauchy = [resultant_forcesN1; resultant_momentN1];
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);
        
    end
    
end

if ~computePostData
    
    %residual (nonlinear)
    rData{1} = + GR2*resultant_forcesN05;
    rData{2} = - V*resultant_forcesN05 + G1*resultant_momentN05;
    rData{3} = - C2*(GammaN1-GammaN)/h + GR2'*(rN1-rN)/h - V'*(phiN1-phiN)/h;
    rData{4} = - C1*(kappaN1-kappaN)/h + G1'*(phiN1-phiN)/h;
    
    %tangent (constant)
    % kData{i,j} = ...
    %
    
end

end


function rot_mat = get_planar_rotation_matrix(angle)

rot_mat = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];

end

function [EA, EI, kGA] = get_material_parameters(materialObject)

if isfield(materialObject,'E') && isfield(materialObject,'A')
    E = materialObject.E;
    area = materialObject.A;
    EA = E * area;
elseif isfield(materialObject,'EA')
    EA = materialObject.EA;
else
    error("No valid input of inertial parameters")
end

if isfield(materialObject,'E') && isfield(materialObject,'I')
    E = materialObject.E;
    inertia = materialObject.I;
    EI = E * inertia;
elseif isfield(materialObject,'EI')
    EI = materialObject.EI;
else
    error("No valid input of inertial parameters")
end

if isfield(materialObject,'shearCorrectionCoefficient') && isfield(materialObject,'G') && isfield(materialObject,'A')
    k = materialObject.shearCorrectionCoefficient;
    area = materialObject.A;
    G = materialObject.G;
    kGA = k * G * area;
elseif isfield(materialObject,'kGA')
    kGA = materialObject.kGA;
else
    error("No valid input of inertial parameters")
end

end