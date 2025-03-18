function [rData, kData, elementEnergy, array] = displacementSCANNDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCANNMOONEYRIVLINDISCRETEGRADIENT Element routine of class solidClass.
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
% In particular algorthmic derivatives for IC, IIC and c are employed. 
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
numberOfDOFs = dimension*size(edof, 2);
I = eye(dimension);
% initialize elementEnergy
elementEnergy.strainEnergy = 0;
% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 3)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;
% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
%     [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    FN = edN * dN_X_I';
    FN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN1 = BMatrix(dN_X_I, FN1);
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
    GN05 = 0.5 * wedge(CN05, CN05); %mid-configuration
    % Determinant
    cN = det(CN);
    cN1 = det(CN1);
    cN05 = det(CN05);
    % invariants for ANN model
    IN = trace(CN);
    IIN = trace(GN);
    IN1 = trace(CN1);
    IIN1 = trace(GN1);    
    % first derivative
    dIC_C = I;
    dIIC_G = I;    
    [~,       ~,        ~,       ~,       ~,           dJ_cN1,  dJstar_cN1] = ANNObject.computeDiffEnergyIIIcANN(CN1,GN1,cN1);
    % second derivative (midpoint)
    d2W_IN05 = 0.5*ANNObject.compute2ndDiffEnergyANN(CN05,GN05,cN05);    
    % discrete gradient
    deltaI = IN1 - IN;
    normDeltaI = abs(deltaI);
    deltaII = IIN1 - IIN;
    normDeltaII = abs(deltaII);
    deltac = cN1 - cN;
    normDeltac = abs(deltac);
    if ~flagNumericalTangent
        discreteGradientThreshold = 1e-10; % FIXME: include into solidObject
        discreteGradientCondition_I = (normDeltaI > discreteGradientThreshold);
        discreteGradientCondition_II = (normDeltaII > discreteGradientThreshold);
        discreteGradientCondition_c = (normDeltac > discreteGradientThreshold);
        if (setupObject.newton.step(setupObject.timeStep) > 1) && (~discreteGradientCondition_I || ~discreteGradientCondition_II || ~discreteGradientCondition_c)
            % disp(1)
        end
    else
        discreteGradientCondition_I = flagDiscreteGradient(k, 1);
        discreteGradientCondition_II = flagDiscreteGradient(k, 2);
        discreteGradientCondition_c = flagDiscreteGradient(k, 3);
    end
    D2W = zeros(4);
    if ~discreteGradientCondition_I || ~discreteGradientCondition_II || ~discreteGradientCondition_c
        [dW_IN05, dW_IIN05, dW_cN05, dW_JN05, dW_JstarN05, dJ_cN05, dJstar_cN05] = ANNObject.computeDiffEnergyIIIcANN(CN05,GN05,cN05);
    end
    if discreteGradientCondition_I
        flagDiscreteGradient(k, 1) = 1;
        tempDW_IN1 = ANNObject.computeEnergyANN(CN1,GN,cN)-ANNObject.computeEnergyANN(CN,GN,cN);
        tempDW_IN = ANNObject.computeEnergyANN(CN1,GN1,cN1)-ANNObject.computeEnergyANN(CN,GN1,cN1);
        DW_I = 1/2*(tempDW_IN1+tempDW_IN)/deltaI;
        tempD2W_IN1a = ANNObject.computeDiffEnergyANNVec(CN1,GN,cN);
        tempD2W_IN1b = ANNObject.computeDiffEnergyANNVec(CN,GN,cN);
        tempD2W_IN1a(2:4) = 0;
        tempD2W_IN1b(:) = 0;
        tempD2W_IN1 = tempD2W_IN1a-tempD2W_IN1b;
        tempD2W_INa = ANNObject.computeDiffEnergyANNVec(CN1,GN1,cN1);
        tempD2W_INb = ANNObject.computeDiffEnergyANNVec(CN,GN1,cN1);
        tempD2W_INb(1) = 0;
        tempD2W_IN = tempD2W_INa-tempD2W_INb;
        D2W(1,:) = 1/2*(tempD2W_IN1+tempD2W_IN)/deltaI;
        D2W(1,1) = D2W(1,1) - 1/2*(tempDW_IN1+tempDW_IN)/(deltaI^2);
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_I = dW_IN05;
        D2W(1,:) = d2W_IN05(1,:);
    end
    if discreteGradientCondition_II
        flagDiscreteGradient(k, 2) = 1;
        tempDW_IIN1 = ANNObject.computeEnergyANN(CN1,GN1,cN)-ANNObject.computeEnergyANN(CN1,GN,cN);
        tempDW_IIN = ANNObject.computeEnergyANN(CN,GN1,cN1)-ANNObject.computeEnergyANN(CN,GN,cN1);
        DW_II = 1/2*(tempDW_IIN1+tempDW_IIN)/deltaII;
        tempD2W_IIN1a = ANNObject.computeDiffEnergyANNVec(CN1,GN1,cN);
        tempD2W_IIN1b = ANNObject.computeDiffEnergyANNVec(CN1,GN,cN);
        tempD2W_IIN1a(3:4) = 0;
        tempD2W_IIN1b(2:4) = 0;
        tempD2W_IIN1 = tempD2W_IIN1a-tempD2W_IIN1b;
        tempD2W_IINa = ANNObject.computeDiffEnergyANNVec(CN,GN1,cN1);
        tempD2W_IINb = ANNObject.computeDiffEnergyANNVec(CN,GN,cN1);
        tempD2W_IINa(1) = 0;
        tempD2W_IINb(1:2) = 0;
        tempD2W_IIN = tempD2W_IINa-tempD2W_IINb;
        D2W(2,:) = 1/2*(tempD2W_IIN1+tempD2W_IIN)/deltaII;
        D2W(2,2) = D2W(2,2) - 1/2*(tempDW_IIN1+tempDW_IIN)/(deltaII^2);        
    else
        flagDiscreteGradient(k, 2) = 0;
        DW_II = dW_IIN05;
        D2W(2,:) = d2W_IN05(2,:);
    end
    if discreteGradientCondition_c
        flagDiscreteGradient(k, 3) = 1;
        tempDW_cN1 = ANNObject.computeEnergyANN(CN1,GN1,cN1)-ANNObject.computeEnergyANN(CN1,GN1,cN);
        tempDW_cN = ANNObject.computeEnergyANN(CN,GN,cN1)-ANNObject.computeEnergyANN(CN,GN,cN);
        DW_c = 1/2*(tempDW_cN1+tempDW_cN)/deltac;
        tempD2W_cN1a = ANNObject.computeDiffEnergyANNVecc(CN1,GN1,cN1);
        tempD2W_cN1b = ANNObject.computeDiffEnergyANNVecc(CN1,GN1,cN);
        tempD2W_cN1b(3) = 0;
        tempD2W_cN1 = tempD2W_cN1a-tempD2W_cN1b;
        tempD2W_cNa = ANNObject.computeDiffEnergyANNVecc(CN,GN,cN1);
        tempD2W_cNb = ANNObject.computeDiffEnergyANNVecc(CN,GN,cN);
        tempD2W_cNa(1:2) = 0;
        tempD2W_cNb(:) = 0;
        tempD2W_cN = tempD2W_cNa-tempD2W_cNb;
        D2W(3,1:3) = 1/2*(tempD2W_cN1+tempD2W_cN)/deltac;
        D2W(3,3) = D2W(3,3) - 1/2*(tempDW_cN1+tempDW_cN)/(deltac^2);        
    else
        flagDiscreteGradient(k, 3) = 0;
        DW_c = dW_cN05;
        D2W(3,:) = d2W_IN05(3,:);        
        D2W(4,:) = d2W_IN05(4,:);        
    end    
    SN05 = 2*(DW_I*dIC_C + wedge(DW_II*dIIC_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1));
    SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];
    if ~computePostData
        %% Residual
        RX = RX + BN05' * SN05_v * detJ * gaussWeight(k);
        %% Tangent
        % geometric tangent
        A1 = 0.5*dN_X_I' * SN05 * dN_X_I;
        geometricTangent = zeros(numberOfDOFs);
        for g = 1:dimension
            geometricTangent(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        % Derivative of dW_C = dW_IC*dIC_C
        KmatCN1 = 0;
        KmatCN05 = 0;
        if ~flagDiscreteGradient(k, 1)
            tempC = wedge(D2W(1,2)*dIIC_G, CN05);
            KmatCN05 = D2W(1,1)*dIC_C([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + dIC_C([1;5;9;4;8;7])*tempC([1;5;9;4;8;7]') + D2W(1,3)*dIC_C([1;5;9;4;8;7])*dJ_cN05*GN05([1;5;9;4;8;7]') + D2W(1,4)*dIC_C([1;5;9;4;8;7])*dJstar_cN05*GN05([1;5;9;4;8;7]');
        else
            tempC = wedge(D2W(1,2)*dIIC_G, CN1);
            KmatCN1 = D2W(1,1)*dIC_C([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + dIC_C([1;5;9;4;8;7])*tempC([1;5;9;4;8;7]') + D2W(1,3)*dIC_C([1;5;9;4;8;7])*dJ_cN1*GN1([1;5;9;4;8;7]') + D2W(1,4)*dIC_C([1;5;9;4;8;7])*dJstar_cN1*GN1([1;5;9;4;8;7]');
        end
        % geometric derivative of wedge(DW_II*dIIC_G,CAlgo)
        KmatGgeoN1 = 0.5*secDiffOperator(DW_II*dIIC_G);
        % material derivative of wedge(dW_IIC*dIIC_G,C)=wedge(C,dW_IIC*dIIC_G)
        temp2 = wedge(CAlgo,dIIC_G);
        KmatGmatN1 = 0;
        KmatGmatN05 = 0;
        if ~flagDiscreteGradient(k, 2)
            temp3N05 = wedge(dIIC_G,CN05);
            KmatGmatN05 = temp2([1;5;9;4;8;7])*D2W(2,1)*dIC_C([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,2)*temp3N05([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,3)*dJ_cN05*GN05([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,4)*dJstar_cN05*GN05([1;5;9;4;8;7]');
        else
            temp3N1 = wedge(dIIC_G,CN1);
            KmatGmatN1 = temp2([1;5;9;4;8;7])*D2W(2,1)*dIC_C([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,2)*temp3N1([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,3)*dJ_cN1*GN1([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*D2W(2,4)*dJstar_cN1*GN1([1;5;9;4;8;7]');
        end
        % geometric derivative of dW_c*G, G =  1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1), GAlgo1 = 0.5 * (GN + GN1); GN1 = 0.5 * wedge(CN1, CN1);  GN = 0.5 * wedge(CN, CN);
        KmatcgeoN1 = DW_c * 2 / 3 *(secDiffOperator(CAlgo) + 0.5*secDiffOperator(CN1));
        % material derivative of dW_c*(1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1)), dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c
        Gtemp = 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1);
        Gtemp_v = [Gtemp(1, 1); Gtemp(2, 2); Gtemp(3, 3); Gtemp(1, 2); Gtemp(2, 3); Gtemp(1, 3)];
        KmatcmatN1 = 0;
        KmatcmatN05 = 0;
        if ~flagDiscreteGradient(k, 3)
            GN05_v = [GN05(1, 1); GN05(2, 2); GN05(3, 3); GN05(1, 2); GN05(2, 3); GN05(1, 3)];
            tempG = wedge(D2W(3,2) * dJ_cN05 * dIIC_G, CN05);
            Kmatcmat1 = 0.5*dW_JN05*Gtemp([1;5;9;4;8;7])*(-1/4*cN05^(-3/2))*GN05([1;5;9;4;8;7]') + 0.5*dW_JstarN05*Gtemp([1;5;9;4;8;7])*(1/4*cN05^(-3/2))*GN05([1;5;9;4;8;7]') + D2W(3,1) * dJ_cN05 * Gtemp([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + Gtemp([1;5;9;4;8;7])*tempG([1;5;9;4;8;7]') + D2W(3,3) * dJ_cN05 * dJ_cN05 * (Gtemp_v * GN05_v') + D2W(3,4) * dJ_cN05 * dJstar_cN05 * (Gtemp_v * GN05_v');
            tempGstar = wedge(D2W(4,2) * dJstar_cN05 * dIIC_G, CN05);
            Kmatcmat2 = D2W(4,1) * dJstar_cN05 * Gtemp([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + Gtemp([1;5;9;4;8;7])*tempGstar([1;5;9;4;8;7]') + D2W(4,3) * dJstar_cN05 * dJ_cN05 * (Gtemp_v * GN05_v') + D2W(4,4) * dJstar_cN05 * dJstar_cN05 * (Gtemp_v * GN05_v');
            KmatcmatN05 = Kmatcmat1 + Kmatcmat2;
        else
            GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3)];
            tempG = wedge(D2W(3,2) * dIIC_G, CN1);
            KmatcmatN1 = D2W(3,1) * Gtemp([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + Gtemp([1;5;9;4;8;7])*tempG([1;5;9;4;8;7]') + D2W(3,3) * (Gtemp_v * GN1_v');
        end
        % Assembly of elasticity tensor
        materialTangentN05 = 4 * (KmatCN05 + KmatGmatN05 + KmatcmatN05);
        materialTangentN1 = 4 * (0.5*KmatcgeoN1 + KmatCN1 + KmatGmatN1 + KmatGgeoN1 + KmatcmatN1);
        KXX = KXX + geometricTangent * detJ * gaussWeight(k);
        KXX = KXX + BN05' * (materialTangentN05 * BN05 + materialTangentN1 * BN1) * detJ * gaussWeight(k);
        % Strain energy
        WN1 = ANNObject.computeEnergyANN(CN1,GN1,cN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + WN1 * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        [dW_IN1, dW_IIN1, dW_cN1] = ANNObject.computeDiffEnergyANN(CN1,GN1,cN1);
        SN1 = 2*(dW_IN1*eye(3) + wedge(dW_IIN1*eye(3), CN1) + dW_cN1 * GN1);
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
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

%% numerical tangent
% if 1
%     %GN05
%     temp = 2;
%     if temp == 0
%         SN052 = 2*(DW_I*dIC_C + wedge(DW_II*dIIC_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1));
%     elseif temp == 1
%         SN052 = 2*(DW_I*dIC_C);
%     elseif temp == 2
%         SN052 = 2*(wedge(DW_II*dIIC_G, CAlgo));
%     elseif temp == 3
%         SN052 = DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1);
%     end
%     %
%     SN05_v2 = [SN052(1, 1); SN052(2, 2); SN052(3, 3); SN052(1, 2); SN052(2, 3); SN052(1, 3)];
%     SN05 = 2*(0*DW_I*dIC_C + 0*wedge(DW_II*dIIC_G, CAlgo) + DW_c * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1));
%     SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];
%     r = BN05'*SN05_v;
%     r2 = BN05'*SN05_v2;
%     h = 1e-6;
%     kNum = zeros(numel(edN1));
%     kNum2 = zeros(numel(edN1));
%     for ii = 1:numel(edN1)
%         edN1Num = edN1;
%         edN1Num(ii) = edN1Num(ii) + h;
%         edN05Num = 0.5 * (edN + edN1Num);
%         FN1Num = edN1Num * dN_X_I';
%         FN05Num = edN05Num * dN_X_I';
%         CN1Num = FN1Num' * FN1Num;
%         CN05Num = FN05Num' * FN05Num;
%         CAlgoNum = 0.5 * (CN + CN1Num);
%         GN1Num = 0.5 * wedge(CN1Num, CN1Num);
%         GN05Num = 0.5 * wedge(CN05Num, CN05Num);
%         GAlgo1Num = 0.5 * (GN + GN1Num); %averaged
%         cN05Num = det(CN05Num);
%         % only midpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         [DW_INum, DW_IINum, DW_cNum] = ANNObject.computeDiffEnergyIIIcANN(CN05Num,GN05Num,cN05Num);
%         % % % % % % % % % % % %
%         SN05Num = 2*(0*DW_INum*dIC_C + 0*wedge(DW_IINum*dIIC_G, CAlgoNum) + DW_cNum * 1 / 3 * (wedge(CAlgoNum, CAlgoNum) + GAlgo1Num));
%         %                 SN05Num = 2*(0*DW_INum*dIC_C + 0*wedge(DW_IINum*dIIC_G, CAlgoNum) + DW_c * 1 / 3 * (wedge(CAlgoNum, CAlgoNum) + GAlgo1Num));
%         %                 SN05Num = SN05Num + 2*DW_cNum * 1 / 3 * (wedge(CAlgo, CAlgo) + GAlgo1);
% 
%         % % % % % % % % % % % %
%         SN05_vNum = [SN05Num(1, 1); SN05Num(2, 2); SN05Num(3, 3); SN05Num(1, 2); SN05Num(2, 3); SN05Num(1, 3)];
%         rNum = BN05' * SN05_vNum;
%         kNum(:,ii) = (rNum-r)/h;
%         % % % % % % % % % % % % % % % % % % % % % % % % % % %
%         if temp == 0
%             SN05Num2 = 2*(DW_INum*dIC_C + wedge(DW_IINum*dIIC_G, CAlgoNum) + DW_cNum * 1 / 3 * (wedge(CAlgoNum, CAlgoNum) + GAlgo1Num));
%         elseif temp == 1
%             SN05Num2 = 2*(DW_INum*dIC_C);
%         elseif temp == 2
%             SN05Num2 = 2*wedge(DW_IINum*dIIC_G, CAlgoNum);
%         elseif temp == 3
%             SN05Num2 = 2*(DW_cNum * 1 / 3 * (wedge(CAlgoNum, CAlgoNum) + GAlgo1Num));
%         end
%         SN05_vNum2 = [SN05Num2(1, 1); SN05Num2(2, 2); SN05Num2(3, 3); SN05Num2(1, 2); SN05Num2(2, 3); SN05Num2(1, 3)];
%         %
%         rNum2 = BN05' * SN05_vNum2;
%         kNum2(:,ii) = (rNum2-r2)/h;
%     end
%     KXX = KXX + geometricTangent * detJ * gaussWeight(k);
%     KXX = KXX + BN05' * (materialTangentN05 * BN05 + materialTangentN1 * BN1) * detJ * gaussWeight(k);
%     %               KXX = KXX + kNum* detJ * gaussWeight(k);
%     %               k2 = BN05' * materialTangentN05temp * BN05 + BN05' * materialTangentN1temp * BN1;
%     %               kNum2-k2
% end
