function [rData, kData, elementData, array] = displacementPhIrreducibleGeometricallyExactHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, flagNumericalTangent)
% Geometrically exact beam in EM conserving, displacement Formulation by
% Stander and Stein

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

% material data
[EA, EI, kGA] = get_material_parameters(materialObject);

% matrices
skew_matrix = [0 , 1;
    -1, 0];
material_matrix = [EA, 0;
    0, kGA];

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% initialize elementEnergy
elementData.strainEnergy = 0;

% residual
rData{2} = zeros(2,1);
rData{1} = zeros(4,1);

%% Gauss loop
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I    = computedN_X_I(dN_xi_k_I, J, k);
    
    % shape function matrices
    dsPhi1 = dN_X_I;
    Phi1   = N_k_I(k, :);
    dsPhi2 = [dN_X_I(1) 0         dN_X_I(2) 0        ;
        0         dN_X_I(1) 0         dN_X_I(2)];
    
    
    GammaN = obj.historyN(e,k).dilatationAndShear;
    
    if ~computePostData
        
        % rotation matrix
        RN05 = get_planar_rotation_matrix(Phi1*phiN05);
        
        % strain quantities
        kappaN05 = dsPhi1*phiN05;
        kappaN1  = dsPhi1*phiN1;
        GammaN1 = GammaN + Phi1*(phiN1-phiN) * RN05' * skew_matrix * (dsPhi2*rN05) + RN05' * dsPhi2*(rN1 - rN);
        GammaN05 = 1/2 * (GammaN + GammaN1);
        
        % update internal variable
        if ~flagNumericalTangent
            obj.historyN1(e,k).dilatationAndShear = GammaN1;
        end
        
        % stress resultants
        resultant_momentN05 = EI*kappaN05;
        resultant_forcesN05 = material_matrix*GammaN05;
        
        % coefficient matrices
        G1  = dsPhi1' * detJ * gaussWeight(k);
        GR2 = dsPhi2' * RN05 * detJ * gaussWeight(k) ;
        V   = Phi1'   * (dsPhi2*rN05)' * skew_matrix * RN05 * detJ * gaussWeight(k);
        
        %residual (nonlinear)
        rData{1} = rData{1} + GR2*resultant_forcesN05;
        rData{2} = rData{2} - V*resultant_forcesN05 + G1*resultant_momentN05;
        
        % potential energy due to strain
        elementData.strainEnergy = elementData.strainEnergy + (1/2 * EI * kappaN1^2 + 1/2 * (GammaN1)'*material_matrix*(GammaN1)) * detJ * gaussWeight(k);
        
    else
        
        RN05 = get_planar_rotation_matrix(Phi1*phiN05);
        GammaN1 = GammaN + Phi1*(phiN1-phiN) * RN05' * skew_matrix * (dsPhi2*rN05) + RN05' * dsPhi2*(rN1 - rN);
        
        % for plotting
        resultant_forcesN1  = material_matrix*GammaN1;
        kappaN1 = dsPhi1*phiN1;
        resultant_momentN1  = EI*kappaN1;
        stressTensor.Cauchy = [resultant_forcesN1; resultant_momentN1];
        array = postStressComputation(array,N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);
        
    end
    
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