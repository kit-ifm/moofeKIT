function [rData, kData, elementEnergy, array] = displacementSCANNEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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

% initialize residual & tangent
RX = rData{1};
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

V = 0;

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    % B-matrix (current configuration)
    BN1 = BMatrix(dN_X_I, FN1);
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1); %cof(CN1)=det(CN1)*CN1^-T
    % Determinant
    cN1 = det(CN1);
    JN1 = sqrt(cN1);
    % first derivative
    dIC_C = I;
    dIIC_G = I;
    [dW_C, dW_G, dW_c, dW_I, dW_II, dW_J, dW_Jstar, dJ_c, dJstar_c] = ANNObject.computeDiffEnergyCGcANN(CN1,GN1,cN1);
%     dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c;
    SN1 = 2*(dW_C + wedge(dW_G, CN1) + dW_c*GN1);
    SN1_v = [SN1(1, 1); SN1(2, 2); SN1(3, 3); SN1(1, 2); SN1(2, 3); SN1(1, 3)];

    % second derivative
    d2W_I_I = ANNObject.compute2ndDiffEnergyANN(CN1,GN1,cN1);
    if ~computePostData
        % Residual
        RX = RX + BN1' * SN1_v * detJ * gaussWeight(k);
        V = V + JN1*detJ*gaussWeight(k);
        if 1
            % Tangent
            % geometric tangent
            A1 = dN_X_I' * SN1 * dN_X_I;
            geometricTangent = zeros(numberOfDOFs);
            for g = 1:dimension
                geometricTangent(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
            end
            % Derivative of dWdev_C = dWdev_IC*dIC_C
            tempC = wedge(d2W_I_I(1,2)*dIIC_G, CN1);
            KmatC = d2W_I_I(1,1)*dIC_C([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + dIC_C([1;5;9;4;8;7])*tempC([1;5;9;4;8;7]') + d2W_I_I(1,3)*dIC_C([1;5;9;4;8;7])*dJ_c*GN1([1;5;9;4;8;7]') + d2W_I_I(1,4)*dIC_C([1;5;9;4;8;7])*dJstar_c*GN1([1;5;9;4;8;7]');
            % geometric derivative of wedge(dW_G,C), dW_G = dW_IIC*dIIC_G
            KmatGgeo = secDiffOperator(dW_G);
            % material derivative of wedge(dW_G,C)=wedge(C,dW_G), dW_G = dW_IIC*dIIC_G
            temp2 = wedge(CN1,dIIC_G);
            temp3 = wedge(dIIC_G,CN1);
            KmatGmat = temp2([1;5;9;4;8;7])*d2W_I_I(2,1)*dIC_C([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_I(2,2)*temp3([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_I(2,3)*dJ_c*GN1([1;5;9;4;8;7]') + temp2([1;5;9;4;8;7])*d2W_I_I(2,4)*dJstar_c*GN1([1;5;9;4;8;7]');
            % geometric derivative of dW_c*G, G = 0.5*wedge(CN1, CN1)
            Kmatcgeo = dW_c * secDiffOperator(CN1);
            % material derivative of dW_c*G, dW_c = dW_J*dJ_c + dW_Jstar*dJstar_c
            GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3)];
            tempG = wedge(d2W_I_I(3,2) * dJ_c * dIIC_G, CN1);
            Kmatcmat1 = dW_J*GN1([1;5;9;4;8;7])*(-1/4*cN1^(-3/2))*GN1([1;5;9;4;8;7]') + dW_Jstar*GN1([1;5;9;4;8;7])*(1/4*cN1^(-3/2))*GN1([1;5;9;4;8;7]') + d2W_I_I(3,1) * dJ_c * GN1([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + GN1([1;5;9;4;8;7])*tempG([1;5;9;4;8;7]') + d2W_I_I(3,3) * dJ_c * dJ_c * (GN1_v * GN1_v') + d2W_I_I(3,4) * dJ_c * dJstar_c * (GN1_v * GN1_v');
            tempGstar = wedge(d2W_I_I(4,2) * dJstar_c * dIIC_G, CN1);
            Kmatcmat2 = d2W_I_I(4,1) * dJstar_c * GN1([1;5;9;4;8;7])*dIC_C([1;5;9;4;8;7]') + GN1([1;5;9;4;8;7])*tempGstar([1;5;9;4;8;7]') + d2W_I_I(4,3) * dJstar_c * dJ_c * (GN1_v * GN1_v') + d2W_I_I(4,4) * dJstar_c * dJstar_c * (GN1_v * GN1_v');
            Kmatcmat = Kmatcmat1 + Kmatcmat2;
            % Assembly of elasticity tensor
            materialTangent = 4 * (KmatC + KmatGgeo + KmatGmat + Kmatcgeo + Kmatcmat);
            KXX = KXX + geometricTangent * detJ * gaussWeight(k) + BN1' * materialTangent * BN1 * detJ * gaussWeight(k);
        end
%         % numerical tangent
%         if 1
%             %
%             if temp == 1
%                  SN12 = 2*(dW_C);
%             elseif temp == 2
%                  SN12 = 2*(wedge(dW_G, CN1));
%             elseif temp == 3
%                 SN12 = dW_c*GN1;
%             end
%             %
%             SN1_v2 = [SN12(1, 1); SN12(2, 2); SN12(3, 3); SN12(1, 2); SN12(2, 3); SN12(1, 3)];
%             SN1 = 2*(wedge(dW_G, CN1) + 0*dW_c*GN1);
%             SN1_v = [SN1(1, 1); SN1(2, 2); SN1(3, 3); SN1(1, 2); SN1(2, 3); SN1(1, 3)];
%             r = BN1'*SN1_v;
%             r2 = BN1'*SN1_v2;
%             h = 1e-6;
%             kNum = zeros(numel(edN1));
%             kNum2 = zeros(numel(edN1));
%             for ii = 1:numel(edN1)
%                 edN1Num = edN1;
%                 edN1Num(ii) = edN1Num(ii) + h;
%                 FN1Num = edN1Num * dN_X_I';
%                 BN1Num = BMatrix(dN_X_I, FN1Num);
%                 CN1Num = FN1Num' * FN1Num;
%                 GN1Num = 0.5 * wedge(CN1Num, CN1Num);
%                 ICNum = trace(CN1Num);
%                 IICNum = trace(GN1Num);
%                 cN1Num = det(CN1Num);
%                 JNum = sqrt(cN1Num);
%                 JstarNum = -JNum;
%                 xNum = [ICNum;IICNum;JNum;JstarNum];
%                 h1Num = w1'*xNum + b1;
%                 dh2dh1Num = diag(activation.derivative(h1Num));
%                 dWvol_JNum = 2*alpha*(JNum-1);
%                 dW_INum = w2'*(dh2dh1Num*w1');
%                 dW_INum(3) = dW_INum(3) + dWvol_JNum;
%                 dW_INum(4) = dW_INum(4);
%                 dW_ICNum = dW_INum(1);
%                 dW_IICNum = dW_INum(2);
%                 dW_JNum = dW_INum(3);
%                 dW_JstarNum = dW_INum(4);
%                 dIC_C = I;
%                 dIIC_G = I;
%                 dJ_cNum = 1/2*(cN1Num)^(-1/2);
%                 dJstar_cNum = -dJ_cNum;
%                 dW_CNum = dW_ICNum*dIC_C;
%                 dW_GNum = dW_IICNum*dIIC_G;
%                 dW_cNum = dW_JNum*dJ_cNum + dW_JstarNum*dJstar_cNum;
%                 %
%                 SN1Num = 2*(0*dW_CNum + wedge(dW_GNum, CN1) + 0*dW_c*GN1);
%                 if temp == 1
%                     SN1Num2 = 2*(dW_CNum);
%                 elseif temp == 2
%                     SN1Num2 = 2*(wedge(dW_GNum, CN1Num));
%                 elseif temp == 3
%                     SN1Num2 = 2*(dW_cNum*GN1Num);
%                 end
%                 %
%                 SN1_vNum = [SN1Num(1, 1); SN1Num(2, 2); SN1Num(3, 3); SN1Num(1, 2); SN1Num(2, 3); SN1Num(1, 3)];
%                 SN1_vNum2 = [SN1Num2(1, 1); SN1Num2(2, 2); SN1Num2(3, 3); SN1Num2(1, 2); SN1Num2(2, 3); SN1Num2(1, 3)];
%                 rNum = BN1' * SN1_vNum;
%                 kNum(:,ii) = (rNum-r)/h;
%                 rNum2 = BN1' * SN1_vNum2;
%                 kNum2(:,ii) = (rNum2-r2)/h;
%             end
%             KXX = KXX + BN1' * materialTangent * BN1 * detJ * gaussWeight(k);% + kNum* detJ * gaussWeight(k);
%             KXX = KXX + geometricTangent * detJ * gaussWeight(k);
%             k2 = BN1' * materialTangent * BN1;
%             kNum2-k2
%             KXX = KXX + BN1' * materialTangent * BN1 * detJ * gaussWeight(k);
%         end
        % Strain energy
        W = ANNObject.computeEnergyANN(CN1,GN1,cN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
%         elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(J)) + c / 2 * (sqrt(J) - 1)^2) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        [~, detJStruct, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
%     disp(V)
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
