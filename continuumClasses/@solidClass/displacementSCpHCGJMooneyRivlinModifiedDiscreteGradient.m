function [rData, kData, elementEnergy, array] = displacementSCpHCGJMooneyRivlinModifiedDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
alpha = materialObject.alpha;
beta = materialObject.beta;
gamma = materialObject.gamma;
epsilon1 = materialObject.epsilon1;
epsilon2 = materialObject.epsilon2;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5*(edN + edN1);

% compute nodal velocities
% edvN1 = obj.vN1(edof, 1:dimension).'; % does not work
DT = setupObject.timeStepSize;
% edvN1 = (edN1-edN)/DT;
edvN05 = (edN1-edN)/DT;

% initialize residual & tangent
RX = rData{1};
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% initialize flagDiscreteGradient
numericalTangentObject = obj.numericalTangentObject;
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 1)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % Deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FxN05);
    % Right Cauchy-Green tensor
    CxN1 = FxN1.' * FxN1;
    CxN05 = FxN05.' * FxN05;
    % Cofactor
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    GxN05 = 0.5 * wedge(CxN05, CxN05);
    % Third invariant
    JxN1 = det(FxN1);
    JxN05 = det(FxN05);
    % mixed/history variables
    FvN05 = edvN05 * dN_X_I';
    CvN05 = FxN05.'*FvN05 + FvN05.'*FxN05;
    GvN05 = wedge(CxN05,CvN05);
    HxN05 = 0.5 * wedge(FxN05, FxN05);
    CN = obj.historyN(e,k).C;
    GN = obj.historyN(e,k).G;
    JN = obj.historyN(e,k).J;  
    CN1 = CN + DT*CvN05;
    GN1 = GN + DT*GvN05;
    JN1 = JN + DT*HxN05(:).'*FvN05(:);
    % update history variables
    if ~flagNumericalTangent
        obj.historyN1(e,k).C = CN1;
        obj.historyN1(e,k).G = GN1;
        obj.historyN1(e,k).J = JN1;
    end
    CN05 = 0.5*(CN + CN1);
    GN05 = 0.5*(GN + GN1);
    JN05 = 0.5*(JN + JN1);
    % Derivative of the strain energy function
    dW_C = alpha * trace(CN1) * I;
    dW_G = beta * trace(GN1) * I;
    dW_J = - gamma*JN1^(-1) + 2*epsilon1*epsilon2*(JN1^(2*epsilon2-1) - JN1^(-2*epsilon2-1));  
    % deltaC = CN1 - CN;
    % normDeltaC = norm(deltaC);
    % if ~flagNumericalTangent
    %     discreteGradientCondition_C = (normDeltaC > 1e-10);
    % else
    %     discreteGradientCondition_C = flagDiscreteGradient(k, 1);
    % end
    % if discreteGradientCondition_C
    %     flagDiscreteGradient(k, 1) = 1;
    %     DW_C = alpha*trace(CN05)*I;
    % else
    %     flagDiscreteGradient(k, 1) = 0;
    %     DW_C = c*(JN05 - 1) - d*JN05^(-1);
    % end
    % deltaG = GN1 - GN;
    % normDeltaG = norm(deltaG);
    % if ~flagNumericalTangent
    %     discreteGradientCondition_G = (normDeltaG > 1e-10);
    % else
    %     discreteGradientCondition_G = flagDiscreteGradient(k, 2);
    % end
    % if discreteGradientCondition_G
    %     flagDiscreteGradient(k, 2) = 1;
    %     DW_G = beta*trace(GN05)*I;
    % else
    %     flagDiscreteGradient(k, 2) = 0;
    %     DW_G = c*(JN05 - 1) - d*JN05^(-1);
    % end
    DW_C = alpha*trace(CN05)*I;
    DW_G = beta*trace(GN05)*I;
    deltaJ = JN1 - JN;
    normDeltaJ = abs(deltaJ);
    if ~flagNumericalTangent
        discreteGradientCondition_J = (normDeltaJ > 1e-6);
    else
        discreteGradientCondition_J = flagDiscreteGradient(k, 1);
    end
    if discreteGradientCondition_J
        flagDiscreteGradient(k, 1) = 1;
        DW_J = (-gamma*(log(JN1)-log(JN))+epsilon1*(JN1^(2*epsilon2)-JN^(2*epsilon2)+JN1^(-2*epsilon2)-JN^(-2*epsilon2)))/deltaJ;
    else
        flagDiscreteGradient(k, 1) = 0;
        % DW_J = c*(JN05 - 1) - d;
        DW_J = - gamma*JN05^(-1) + 2*epsilon1*epsilon2*(JN05^(2*epsilon2-1) - JN05^(-2*epsilon2-1));
    end
    % Second Piola Kirchhoff stress tensor
    SN1 = 2 * (dW_C + wedge(dW_G, CxN1) + 1/2*dW_J*JxN1^(-1)*GxN1);
    SAlgo = 2 * (DW_C + wedge(DW_G, CxN05) + 1/2*DW_J*JxN05^(-1)*GxN05);
    SAlgo_v = matrixToVoigt(SAlgo, 'stress');
    if ~computePostData
        % Residual
        RX = RX + BN05.' * SAlgo_v * detJ * gaussWeight(k);
        % Strain energy
        W = alpha/2 * ((trace(CN1))^2-9) + beta/2 * ((trace(GN1))^2 - 9) - gamma*log(JN1) + epsilon1*(JN1^(2*epsilon2)+JN1^(-2*epsilon2)-2);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
        % % Tangent
        % Sigma_I_I = c/2*JxN1^(-2) + d*JxN1^(-3);    
        % % Derivative of wedge(Sigma_G,CN05)
        % if flagSymmetric
        %     Kmat1 = secDiffOperator(Sigma_G);
        % else
        %     Kmat1 = b*[  0     1     1     0     0     0     0     0     0;
        %                  1     0     1     0     0     0     0     0     0;
        %                  1     1     0     0     0     0     0     0     0;
        %                  0     0     0     0     0     0    -1     0     0;
        %                  0     0     0     0     0     0     0    -1     0;
        %                  0     0     0     0     0     0     0     0    -1;
        %                  0     0     0    -1     0     0     0     0     0;
        %                  0     0     0     0    -1     0     0     0     0;
        %                  0     0     0     0     0    -1     0     0     0];
        % end
        % % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I
        % if flagSymmetric
        %     Kmat2 = Sigma_I * secDiffOperator(CN1);
        % else
        %     Kmat2 = Sigma_I*[   0     CN1(3,3)     CN1(2,2)     0     -CN1(3,2)     0     0     -CN1(2,3)     0;
        %                          CN1(3,3)     0     CN1(1,1)     0     0     -CN1(3,1)     0     0     -CN1(1,3);
        %                          CN1(2,2)     CN1(1,1)     0     -CN1(2,1)     0     0     -CN1(1,2)     0     0;
        %                          0     0     -CN1(2,1)     0     CN1(3,1)     0    -CN1(3,3)     0     CN1(2,3);
        %                          -CN1(3,2)     0     0     CN1(3,1)     0     0     0    -CN1(1,1)     CN1(1,2);
        %                          0     -CN1(3,1)     0     0     0     0     CN1(3,2)     CN1(2,1)    -CN1(2,2);
        %                          0     0     -CN1(1,2)    -CN1(3,3)     0     CN1(3,2)     0     CN1(1,3)     0;
        %                          -CN1(2,3)     0     0     0    -CN1(1,1)     CN1(2,1)     CN1(1,3)     0     0;
        %                          0     -CN1(1,3)     0     CN1(2,3)     CN1(1,2)    -CN1(2,2)     0     0     0];
        % end
        % % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
        % if flagSymmetric
        %     GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3)];
        % else
        %     GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3); GN1(2, 1); GN1(3, 2); GN1(3, 1)];
        % end
        % Kmat3 = Sigma_I_I * 1/2 * JN1^(-1) * (GN1_v * GN1_v');
        % % Assembly of elasticity tensor
        % ELA = 4 * (Kmat1 + Kmat2 + Kmat3);
        % MAT = BN1' * ELA * BN1 * detJ * gaussWeight(k);
        % GEO = kron(dN_X_I' * SN1 * dN_X_I * detJ * gaussWeight(k),eye(3));
        % KXX = KXX + MAT + GEO;
    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;

    % pass flagDiscreteGradient
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
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