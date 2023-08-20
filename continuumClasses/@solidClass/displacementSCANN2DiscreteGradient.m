function [rData, kData, elementEnergy, array] = displacementSCANN2DiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCANNMOONEYRIVLIN2DISCRETEGRADIENT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a trained artificial neural network (ANN) 
% constitutive model which was trained with an isotropic, hyperelastic 
% Mooney-Rivlin ('MooneyRivlin') model.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for dynamic simulations where an energy-momentum 
% consistent integration scheme is used which is based on the conecpt of 
% discrete gradients in the sense of Gonzalez ('DiscreteGradient').
% In particular algorthmic derivatives for C, G and c are employed. 
%
% CALL
% displacementSCANNMooneyRivlinDiscreteGradient(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2018.01.013
%
% SEE ALSO
% displacementSCANNMooneyRivlinEndpoint
% displacementSCANNMooneyRivlinMidpoint
% displacementSCMooneyRivlinMidpoint,
% displacementSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;
% load trained ann model
ANNObject = obj.artificialNeuralNetworkObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
edof = meshObject.edof(e, :);
dimension = obj.dimension;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5 * (edN + edN1);

% initialize residual & tangent
RX = rData{1};
KXX = kData{1, 1};

I = eye(3);

% initialize elementEnergy
elementEnergy.strainEnergy = 0;
% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 3)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;
% Run through all Gauss points
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    FN = edN * dN_X_I';
    FN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FN05);
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    CN = FN' * FN;
    CN05 = FN05' * FN05;
    CAlgo = 0.5 * (CN + CN1); %averaged
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1); %cof(CN1)=det(CN1)*CN1^-T
    GN = 0.5 * wedge(CN, CN);
    GAlgo1 = 0.5 * (GN + GN1); %averaged
    GAlgo2 = 0.5 * wedge(CAlgo, CAlgo); %averaged
    GN05 = 0.5 * wedge(CN05, CN05); %mid-configuration
    % Determinant
    cN = det(CN);
    cN1 = det(CN1);
    cN05 = det(CN05);

    % first derivative
    dIC_C = I;
    dIIC_G = I;
    [dW_CN05, dW_GN05, dW_cN05, dW_IN05, dW_IIN05, dW_JN05, dW_JstarN05, dJ_cN05, dJstar_cN05] = ANNObject.computeDiffEnergyCGcANN(CN05,GN05,cN05);
%     dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c;

    deltaC = CN1 - CN;
    normDeltaC = sqrt(innerProduct(deltaC, deltaC));
    deltaG = GN1 - GN;
    normDeltaG = sqrt(innerProduct(deltaG, deltaG));
    deltac = cN1 - cN;
    normDeltac = abs(deltac);
    
    % discrete gradient
    if ~flagNumericalTangent
        discreteGradientThreshold = 1e-16; % FIXME: include into solidObject
        discreteGradientCondition_C = (normDeltaC > discreteGradientThreshold);
        discreteGradientCondition_G = (normDeltaG > discreteGradientThreshold);
        discreteGradientCondition_c = (normDeltac > discreteGradientThreshold);
        if (setupObject.newton.step > 1) && (~discreteGradientCondition_C || ~discreteGradientCondition_G || ~discreteGradientCondition_c)
            disp(1)
        end
    else
        discreteGradientCondition_C = flagDiscreteGradient(k, 1);
        discreteGradientCondition_G = flagDiscreteGradient(k, 2);
        discreteGradientCondition_c = flagDiscreteGradient(k, 3);
    end
    if discreteGradientCondition_C
        flagDiscreteGradient(k, 1) = 1;
        [dW_C, ~, ~] = ANNObject.computeDiffEnergyCGcANN(CN05,GN,cN);
        dW_CN1 = dW_C + ((ANNObject.computeEnergyANN(CN1,GN,cN)-ANNObject.computeEnergyANN(CN,GN,cN)-sum(sum(dW_C.*deltaC)))/(innerProduct(deltaC, deltaC)))*deltaC;
        [dW_C, ~, ~] = ANNObject.computeDiffEnergyCGcANN(CN05,GN1,cN1);
        dW_CN = dW_C + ((ANNObject.computeEnergyANN(CN1,GN1,cN1)-ANNObject.computeEnergyANN(CN,GN1,cN1)-sum(sum(dW_C.*deltaC)))/(innerProduct(deltaC, deltaC)))*deltaC;
        DW_C = 1/2*(dW_CN1+dW_CN);
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_C = dW_CN05;
    end
    if discreteGradientCondition_G
        flagDiscreteGradient(k, 2) = 1;
        [~, dW_G, ~] = ANNObject.computeDiffEnergyCGcANN(CN1,GN05,cN);
        dW_GN1 = dW_G + ((ANNObject.computeEnergyANN(CN1,GN1,cN)-ANNObject.computeEnergyANN(CN1,GN,cN)-sum(sum(dW_G.*deltaG)))/(innerProduct(deltaG, deltaG)))*deltaG;
        [~, dW_G, ~] = ANNObject.computeDiffEnergyCGcANN(CN,GN05,cN1);
        dW_GN = dW_G + ((ANNObject.computeEnergyANN(CN,GN1,cN1)-ANNObject.computeEnergyANN(CN,GN,cN1)-sum(sum(dW_G.*deltaG)))/(innerProduct(deltaG, deltaG)))*deltaG;
        DW_G = 1/2*(dW_GN1+dW_GN);
    else
        flagDiscreteGradient(k, 2) = 0;
        DW_G = dW_GN05;
    end    
    if discreteGradientCondition_c
        flagDiscreteGradient(k, 3) = 1;
        dW_cN1 = (ANNObject.computeEnergyANN(CN1,GN1,cN1)-ANNObject.computeEnergyANN(CN1,GN1,cN))/deltac;
        dW_cN = (ANNObject.computeEnergyANN(CN,GN,cN1)-ANNObject.computeEnergyANN(CN,GN,cN))/deltac;
        DW_c = 1/2*(dW_cN1+dW_cN);
    else
        flagDiscreteGradient(k, 3) = 0;
        DW_c = dW_cN05;
    end    
%     SN05 = 2*(DW_C + wedge(DW_G, CN05) + DW_c*GN05);
    SN05 = 2*(DW_C + wedge(DW_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1));
%     SN05 = 2*(DW_C + wedge(DW_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo2));
    SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];
    if ~computePostData
        % Residual
        RX = RX + BN05' * SN05_v * detJ * gaussWeight(k);
        % Strain energy
        WN1 = ANNObject.computeEnergyANN(CN1,GN1,cN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + WN1 * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -0.5 * (A(2, 3) + A(3, 2));
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -0.5 * (A(1, 3) + A(3, 1));
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 3) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 4) = -0.5 * A(3, 3);
D(4, 5) = 0.25 * (A(1, 3) + A(3, 1));
D(4, 6) = 0.25 * (A(2, 3) + A(3, 2));
D(5, 1) = -0.5 * (A(2, 3) + A(3, 2));
D(5, 4) = 0.25 * (A(1, 3) + A(3, 1));
D(5, 5) = -0.5 * A(1, 1);
D(5, 6) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 2) = -0.5 * (A(1, 3) + A(3, 1));
D(6, 4) = 0.25 * (A(2, 3) + A(3, 2));
D(6, 5) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 6) = -0.5 * A(2, 2);
end