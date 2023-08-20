function [rData, kData, elementEnergy, array] = displacementSCANNMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCMOONEYRIVLINENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic Mooney-Rivlin
% ('MooneyRivlin') model (attention: nonlinear geometric and material/
% stress-strain relation).
% The formulation is based on a polyconvexity inspired strain-energy
% density function.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementSCMooneyRivlinEndpoint(obj,setupObject,computePostData)
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
% displacementSCMooneyRivlinMidpoint,
% displacementSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;
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

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);


% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    FN = edN * dN_X_I';
    FN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FN05);
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    CN = FN' * FN;
    CN05 = FN05' * FN05; %averaged
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1); %cof(CN1)=det(CN1)*CN1^-T
    GN05 = 0.5 * wedge(CN05, CN05); %mid-configuration
    % Determinant
    cN1 = det(CN1);
    cN05 = det(CN05);
    % first derivative
    dIC_C = I;
    dIIC_G = I;
    [dW_CN05, dW_GN05, dW_cN05, dW_IN05, dW_IIN05, dW_JN05, dW_JstarN05, dJ_cN05, dJstar_cN05] = ANNObject.computeDiffEnergyCGcANN(CN05,GN05,cN05);
%     dW_cN05 = dW_JN05*dJ_cN05 + dW_JstarN05*dJstar_cN05;

    SN05 = 2*(dW_CN05 + wedge(dW_GN05, CN05) + dW_cN05*GN05);
    SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];

    % second derivative
    d2W_I_IN05 = 0.5*ANNObject.compute2ndDiffEnergyANN(CN05,GN05,cN05);
    if ~computePostData
        % Residual
        RX = RX + BN05' * SN05_v * detJ * gaussWeight(k);
        % Tangent
        % geometric tangent
        A1 = 0.5*dN_X_I' * SN05 * dN_X_I;
        geometricTangent = zeros(numberOfDOFs);
        for g = 1:dimension
            geometricTangent(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        % Derivative of dWdev_C = dWdev_IC*dIC_C
        tempC = wedge(d2W_I_IN05(1,2)*dIIC_G, CN05);
        KmatC = d2W_I_IN05(1,1)*dIC_C([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + dIC_C([1;5;9;4;8;7])*tempC([1;5;9;4;8;7]') + d2W_I_IN05(1,3)*dIC_C([1;5;9;4;8;7])*dJ_cN05*GN05([1;5;9;4;8;7]') + d2W_I_IN05(1,4)*dIC_C([1;5;9;4;8;7])*dJstar_cN05*GN05([1;5;9;4;8;7]');
        % geometric derivative of wedge(dW_G,C), dW_G = dW_IIC*dIIC_G
        KmatGgeo = 0.5*secDiffOperator(dW_GN05);
        % material derivative of wedge(dW_G,C)=wedge(C,dW_G), dW_G = dW_IIC*dIIC_G
        temp2 = wedge(CN05,dIIC_G);
        temp3 = wedge(dIIC_G,CN05);
        KmatGmat = temp2([1;5;9;4;8;7])*d2W_I_IN05(2,1)*dIC_C([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_IN05(2,2)*temp3([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_IN05(2,3)*dJ_cN05*GN05([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_IN05(2,4)*dJstar_cN05*GN05([1;5;9;4;8;7]');
        % geometric derivative of dW_c*G, G = 0.5*wedge(CN1, CN1)
        Kmatcgeo = 0.5 * dW_cN05 * secDiffOperator(CN05);
        % material derivative of dW_c*G, dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c
        GN05_v = [GN05(1, 1); GN05(2, 2); GN05(3, 3); GN05(1, 2); GN05(2, 3); GN05(1, 3)];
        tempG = wedge(d2W_I_IN05(3,2) * dJ_cN05 * dIIC_G, CN05);
        Kmatcmat1 = 0.5*dW_JN05*GN05([1;5;9;4;8;7])*(-1/4*cN05^(-3/2))*GN05([1;5;9;4;8;7]') + 0.5*dW_JstarN05*GN05([1;5;9;4;8;7])*(1/4*cN05^(-3/2))*GN05([1;5;9;4;8;7]') + d2W_I_IN05(3,1) * dJ_cN05 * GN05([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + GN05([1;5;9;4;8;7])*tempG([1;5;9;4;8;7]') + d2W_I_IN05(3,3) * dJ_cN05 * dJ_cN05 * (GN05_v * GN05_v') + d2W_I_IN05(3,4) * dJ_cN05 * dJstar_cN05 * (GN05_v * GN05_v');
        tempGstar = wedge(d2W_I_IN05(4,2) * dJstar_cN05 * dIIC_G, CN05);
        Kmatcmat2 = d2W_I_IN05(4,1) * dJstar_cN05 * GN05([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + GN05([1;5;9;4;8;7])*tempGstar([1;5;9;4;8;7]') + d2W_I_IN05(4,3) * dJstar_cN05 * dJ_cN05 * (GN05_v * GN05_v') + d2W_I_IN05(4,4) * dJstar_cN05 * dJstar_cN05 * (GN05_v * GN05_v');
        Kmatcmat = Kmatcmat1 + Kmatcmat2;
        % Assembly of elasticity tensor
        materialTangent = 4 * (KmatC + KmatGgeo + KmatGmat + Kmatcgeo + Kmatcmat);
        KXX = KXX + geometricTangent * detJ * gaussWeight(k);
        KXX = KXX + BN05' * materialTangent * BN05 * detJ * gaussWeight(k);
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

%         % numerical tangent
%         if 0
%             temp = 1;
%             %
%             if temp == 1
%                  SN052 = 2*(dW_CN05);
%             elseif temp == 2
%                  SN052 = 2*(wedge(dW_GN05, CN05));
%             elseif temp == 3
%                 SN052 = dW_cN05*GN05;
%             end
%             %
%             SN05_v2 = [SN052(1, 1); SN052(2, 2); SN052(3, 3); SN052(1, 2); SN052(2, 3); SN052(1, 3)];
%             SN05 = 2*(0*dW_CN05 + 0*wedge(dW_GN05, CN05) + dW_cN05*GN05);
%             SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];
%             r = BN05'*SN05_v;
%             r2 = BN05'*SN05_v2;
%             h = 1e-6;
%             kNum = zeros(numel(edN1));
%             kNum2 = zeros(numel(edN1));
%             for ii = 1:numel(edN1)
%                 edN1Num = edN1;
%                 edN1Num(ii) = edN1Num(ii) + h;
%                 edN05Num = 0.5 * (edN + edN1Num);
%                 FN05Num = edN05Num * dN_X_I';
%                 BN05Num = BMatrix(dN_X_I, FN05Num);
%                 CN05Num = FN05Num' * FN05Num;
%                 GN05Num = 0.5 * wedge(CN05Num, CN05Num);
%                 ICN05Num = trace(CN05Num);
%                 IICN05Num = trace(GN05Num);
%                 cN05Num = det(CN05Num);
%                 JN05Num = sqrt(cN05Num);
%                 JstarN05Num = -JN05Num;
%                 xN05Num = [ICN05Num;IICN05Num;JN05Num;JstarN05Num];
%                 h1N05Num = w1'*xN05Num + b1;
%                 dh2dh1N05Num = diag(activation.derivative(h1N05Num));
%                 dWvol_JNumN05 = 2*alpha*(JN05Num-1);
%                 dW_INumN05 = w2'*(dh2dh1N05Num*w1');
%                 dW_INumN05(3) = dW_INumN05(3) + dWvol_JNumN05;
%                 dW_INumN05(4) = dW_INumN05(4);
%                 dW_ICN05Num = dW_INumN05(1);
%                 dW_IICN05Num = dW_INumN05(2);
%                 dW_JN05Num = dW_INumN05(3);
%                 dW_JstarN05Num = dW_INumN05(4);
%                 dIC_CN05 = I;
%                 dIIC_GN05 = I;
%                 dJ_cN05Num = 1/2*(cN05Num)^(-1/2);
%                 dJstarN05_cNum = -dJ_cN05Num;
%                 dW_CN05Num = dW_ICN05Num*dIC_CN05;
%                 dW_GN05Num = dW_IICN05Num*dIIC_GN05;
%                 dW_cN05Num = dW_JN05Num*dJ_cN05Num + dW_JstarN05Num*dJstarN05_cNum;
%                 %
%                 SN05Num = 2*(0*dW_CN05Num + 0*wedge(dW_GN05Num, CN05Num) + dW_cN05Num*GN05Num);
%                 if temp == 1
%                     SN05Num2 = 2*(dW_CN05Num);
%                 elseif temp == 2
%                     SN05Num2 = 2*(wedge(dW_GN05Num, CN05Num));
%                 elseif temp == 3
%                     SN05Num2 = 2*(dW_cN05Num*GN05Num);
%                 end
%                 %
%                 SN05_vNum = [SN05Num(1, 1); SN05Num(2, 2); SN05Num(3, 3); SN05Num(1, 2); SN05Num(2, 3); SN05Num(1, 3)];
%                 SN05_vNum2 = [SN05Num2(1, 1); SN05Num2(2, 2); SN05Num2(3, 3); SN05Num2(1, 2); SN05Num2(2, 3); SN05Num2(1, 3)];
%                 rNum = BN05' * SN05_vNum;
%                 kNum(:,ii) = (rNum-r)/h;
%                 rNum2 = BN05' * SN05_vNum2;
%                 kNum2(:,ii) = (rNum2-r2)/h;
%             end
%             KXX = KXX + kNum* detJ * gaussWeight(k);
% %             k2 = BN1' * materialTangent * BN1;
% %             kNum2-k2
%         end

