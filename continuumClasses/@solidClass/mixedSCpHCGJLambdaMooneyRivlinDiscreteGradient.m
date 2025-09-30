function [rData, kData, elementEnergy, array] = mixedSCpHCGJLambdaMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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

DT = setupObject.timeStepSize;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);
numberOfXDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalNodes = size(M_k_I, 2);

dimension = obj.dimension;

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5*(edN + edN1);

% compute nodal velocities
% edvN1 = obj.vN1(edof, 1:dimension).'; % does not work
edvN05 = (edN1-edN)/DT;

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

edAlphaN = mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedCNv = edAlphaN(1:6*numberOfInternalNodes);
extractedGNv = edAlphaN(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN = edAlphaN(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes);
extractedGN1v = edAlphaN1(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedLambdaCN1v = edAlphaN1(13*numberOfInternalNodes+1:19*numberOfInternalNodes);
extractedLambdaGN1v = edAlphaN1(19*numberOfInternalNodes+1:25*numberOfInternalNodes);
extractedLambdaJN1 = edAlphaN1(25*numberOfInternalNodes+1:26*numberOfInternalNodes);

% initialize residual
RX = rData{1};
RDotCv = rData{2};
RDotGv = rData{3};
RDotJ = rData{4};
RLambdaCv = rData{5};
RLambdaGv = rData{6};
RLambdaJ = rData{7};

% initialize tangent
KXX = kData{1, 1};
KXLambdaC = kData{1, 5};
KCC = kData{2, 2};
KCG = kData{2, 3};
KCLambdaC = kData{2, 5};
KCLambdaG = kData{2, 6};
KCLambdaJ = kData{2, 7};
KGC = kData{3, 2};
KGLambdaG = kData{3, 6};
KGLambdaJ = kData{3, 7};
KJJ = kData{4, 4};
KJLambdaJ = kData{4, 7};
KLambdaCX = kData{5, 1};
KLambdaCC = kData{5, 2};
KLambdaGC = kData{6, 2};
KLambdaGG = kData{6, 3};
KLambdaJX = kData{7, 1};
KLambdaJC = kData{7, 2};
KLambdaJG = kData{7, 3};
KLambdaJJ = kData{7, 4};

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
    % compute the values of the variables at the current Gauss point
    CNv = reshape(extractedCNv, 6, []) * M_k_I(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    CN1v = reshape(extractedCN1v, 6, []) * M_k_I(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    CN05 = 0.5 * (CN + CN1);

    GNv = reshape(extractedGNv, 6, []) * M_k_I(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M_k_I(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    GN05 = 0.5 * (GN + GN1);

    JN = extractedJN.' * M_k_I(k, :)';
    JN1 = extractedJN1.' * M_k_I(k, :)';
    JN05 = 0.5 * (JN + JN1);

    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M_k_I(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');
    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M_k_I(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');
    lambdaJN1 = extractedLambdaJN1.' * M_k_I(k, :)';
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
    % Determinant
    JxN1 = det(FxN1);
    JxN05 = det(FxN05);
    % Derivative of the strain energy function
    DW_C = a * I;
    DW_G = b * I;
    deltaJ = JN1 - JN;
    normDeltaJ = abs(deltaJ);
    if ~flagNumericalTangent
        discreteGradientCondition_J = (normDeltaJ > 1e-10);
    else
        discreteGradientCondition_J = flagDiscreteGradient(k, 1);
    end
    if discreteGradientCondition_J
        flagDiscreteGradient(k, 1) = 1;
        DW_J = (c/2*(JN1-1)^2 - d*log(JN1) - (c/2*(JN-1)^2 - d*log(JN)))/deltaJ;
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_J = c*(JN05 - 1) - d*JN05^(-1);
    end

    dW_J = c*(JN1 - 1) - d*JN1^(-1);
    % d_J_d_J_W = c + d*JN1^(-2);
    % Second Piola Kirchhoff stress tensor
    SN1 = 2 * (DW_C + wedge(DW_G, CxN1) + 1/2*dW_J*JxN1^(-1)*GxN1);
    SN05 = 2 * (lambdaCN1 + wedge(lambdaGN1, CxN05) + 1/2*lambdaJN1*JxN05^(-1)*GxN05);
    FvN05 = edvN05 * dN_X_I';
    CvN05 = FxN05.'*FvN05 + FvN05.'*FxN05;
    GvN05 = wedge(CxN05,CvN05);
    HxN05 = 0.5 * wedge(FxN05, FxN05);
    SN05_v = matrixToVoigt(SN05, 'stress');
    if ~computePostData
        % Residual
        RX = RX + BN05.' * SN05_v * detJ * gaussWeight(k);
        RDotCv = RDotCv + kron(M_k_I(k, :)',matrixToVoigt((CN1 - CN)/DT - CvN05,'stress')) * detJ * gaussWeight(k);
        RDotGv = RDotGv + kron(M_k_I(k, :)',matrixToVoigt((GN1 - GN)/DT - GvN05,'stress')) * detJ * gaussWeight(k);
        RDotJ = RDotJ + M_k_I(k, :)' * ((JN1 - JN)/DT - HxN05(:).'*FvN05(:)) * detJ * gaussWeight(k);
        RLambdaCv = RLambdaCv + kron(M_k_I(k, :)', matrixToVoigt(DW_C,'stress')-lambdaCN1v) * detJ * gaussWeight(k);
        RLambdaGv = RLambdaGv + kron(M_k_I(k, :)', matrixToVoigt(DW_G, 'stress')-lambdaGN1v) * detJ * gaussWeight(k);
        RLambdaJ = RLambdaJ + M_k_I(k, :)' * (DW_J-lambdaJN1) * detJ * gaussWeight(k);
        % Strain energy
        W = a * (trace(CN1) - 3) + b * (trace(GN1) - 3) + c / 2 * (JN1 - 1)^2 - d * log(JN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RDotCv;
    rData{3} = RDotGv;
    rData{4} = RDotJ;
    rData{5} = RLambdaCv;
    rData{6} = RLambdaGv;
    rData{7} = RLambdaJ;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 5} = KXLambdaC;
    kData{2, 2} = KCC;
    kData{2, 3} = KCG;
    kData{2, 5} = KCLambdaC;
    kData{2, 6} = KCLambdaG;
    kData{2, 7} = KCLambdaJ;
    kData{3, 2} = KGC;
    kData{3, 6} = KGLambdaG;
    kData{3, 7} = KGLambdaJ;
    kData{4, 4} = KJJ;
    kData{4, 7} = KJLambdaJ;
    kData{5, 1} = KLambdaCX;
    kData{5, 2} = KLambdaCC;
    kData{6, 2} = KLambdaGC;
    kData{6, 3} = KLambdaGG;
    kData{7, 1} = KLambdaJX;
    kData{7, 2} = KLambdaJC;
    kData{7, 3} = KLambdaJG;
    kData{7, 4} = KLambdaJJ;

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
