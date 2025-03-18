function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDSCMOONEYRIVLINDISCRETEGRADIENT Element routine of class solidClass.
%
% FORMULATION
% This is a 'mixed'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic Mooney-Rivlin
% ('MooneyRivlin') model (attention: nonlinear geometric and material/
% stress-strain relation).
% The formulation is based on a polyconvexity inspired strain-energy
% density function.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for dynamic simulation where for the concept of
% the energy and momentum consistent discrete gradient in the sense of
% Gonzalez is is used ('DiscreteGradient').
%
% CALL
% mixedSCMooneyRivlinEndpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, C, G, c, LambdaC, LambdaG, Lambdac)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, C, G, c, LambdaC,
%        LambdaG, Lambdac)
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
% https://doi.org/10.1016/j.cma.2018.01.013
%
% SEE ALSO
% mixedSCMooneyRivlinMidpoint
% mixedSCMooneyRivlinEndpoint
%
% CREATOR(S)
% Felix Zaehringer, Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
if strcmpi(mixedFEObject.typeShapeFunction, 'sameOrder')
    M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
    M_C = M_k_I;
    M_G = M_k_I;
    M_c = M_k_I;
    M_LambdaC = M_k_I;
    M_LambdaG = M_k_I;
    M_Lambdac = M_k_I;
elseif strcmpi(mixedFEObject.typeShapeFunction, 'detailedOrder')
    M_C = mixedFEObject.shapeFunctionObject.N{1, 1};
    M_G = mixedFEObject.shapeFunctionObject.N{2, 1};
    M_c = mixedFEObject.shapeFunctionObject.N{3, 1};
    M_LambdaC = mixedFEObject.shapeFunctionObject.N{4, 1};
    M_LambdaG = mixedFEObject.shapeFunctionObject.N{5, 1};
    M_Lambdac = mixedFEObject.shapeFunctionObject.N{6, 1};
end

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;
numberOfXDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalNodes_C = size(M_C, 2);
numberOfInternalNodes_G = size(M_G, 2);
numberOfInternalNodes_c = size(M_c, 2);
numberOfInternalNodes_LambdaC = size(M_LambdaC, 2);
numberOfInternalNodes_LambdaG = size(M_LambdaG, 2);
numberOfInternalNodes_Lambdac = size(M_Lambdac, 2);

dimension = obj.dimension;

I6 = [eye(3), zeros(3); zeros(3), 2 * eye(3)];

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof(e, :), 1:dimension).';

edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
edN05 = 1 / 2 * (edN + edN1);

edAlphaN = mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedCNv = edAlphaN(1:6*numberOfInternalNodes_C);
extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes_C);
extractedGNv = edAlphaN(6*numberOfInternalNodes_C+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G);
extractedGN1v = edAlphaN1(6*numberOfInternalNodes_C+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G);
extractedcN = edAlphaN(6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c);
extractedcN1 = edAlphaN1(6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c);
extractedLambdaCN1v = edAlphaN1(6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+6*numberOfInternalNodes_LambdaC);
extractedLambdaGN1v = edAlphaN1(6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+6*numberOfInternalNodes_LambdaC+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+6*numberOfInternalNodes_LambdaC+6*numberOfInternalNodes_LambdaG);
extractedLambdacN1 = edAlphaN1(6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+6*numberOfInternalNodes_LambdaC+6*numberOfInternalNodes_LambdaG+1:6*numberOfInternalNodes_C+6*numberOfInternalNodes_G+numberOfInternalNodes_c+6*numberOfInternalNodes_LambdaC+6*numberOfInternalNodes_LambdaG+numberOfInternalNodes_Lambdac);

% initialize residual
RX = rData{1};
RCv = rData{2};
RGv = rData{3};
Rc = rData{4};
RLambdaCv = rData{5};
RLambdaGv = rData{6};
RLambdac = rData{7};

% initialize tangent
KXX = kData{1, 1};
KXLambdaC = kData{1, 5};
KCC = kData{2, 2};
KCG = kData{2, 3};
KCLambdaC = kData{2, 5};
KCLambdaG = kData{2, 6};
KCLambdac = kData{2, 7};
KGC = kData{3, 2};
KGLambdaG = kData{3, 6};
KGLambdac = kData{3, 7};
Kcc = kData{4, 4};
KcLambdac = kData{4, 7};
KLambdaCX = kData{5, 1};
KLambdaCC = kData{5, 2};
KLambdaGC = kData{6, 2};
KLambdaGG = kData{6, 3};
KLambdacX = kData{7, 1};
KLambdacC = kData{7, 2};
KLambdacG = kData{7, 3};
KLambdacc = kData{7, 4};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, numberOfGausspoints, 1)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute the values of the variables at the current Gauss point
    CNv = reshape(extractedCNv, 6, []) * M_C(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    CN1v = reshape(extractedCN1v, 6, []) * M_C(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    CN05 = 0.5 * (CN + CN1);

    GNv = reshape(extractedGNv, 6, []) * M_G(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M_G(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    GN05 = 0.5 * (GN + GN1);

    cN = extractedcN.' * M_c(k, :)';
    cN1 = extractedcN1.' * M_c(k, :)';
    cN05 = 0.5 * (cN + cN1);

    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M_LambdaC(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');

    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M_LambdaG(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');

    lambdacN1 = extractedLambdacN1.' * M_Lambdac(k, :)';

    % deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN05 = edN05 * dN_X_I';

    % strain measures
    CxN1 = FxN1.' * FxN1;
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    cxN1 = det(CxN1);

    % nodal operator matrices
    BN1 = BMatrix(dN_X_I, FxN1);
    BN05 = BMatrix(dN_X_I, FxN05);

    if ~computePostData
        % ENERGY
        % elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3) - d * log(sqrt(cxN1)) + c / 2 * (sqrt(cxN1) - 1)^2) * detJ * gaussWeight(k);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2) * detJ * gaussWeight(k);

        % DERIVATIVES
        DW_C = a * eye(3);
        DW_G = b * eye(3);
        deltac = cN1 - cN;
        normDeltac = abs(deltac);
        if ~flagNumericalTangent
            discreteGradientCondition_c = (normDeltac > 1e-10);
        else
            discreteGradientCondition_c = flagDiscreteGradient(k, 1);
        end
        if discreteGradientCondition_c
            flagDiscreteGradient(k, 1) = 1;
            DW_c = (-d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2 + d * log(sqrt(cN)) - c / 2 * (sqrt(cN) - 1)^2) / deltac;
        else
            flagDiscreteGradient(k, 1) = 0;
            DW_c = -d / (2 * cN05) + c / 2 * (1 - 1 / sqrt(cN05));
        end

        % RESIDUAL
        RX = RX + 2 * BN05.' * lambdaCN1v * detJ * gaussWeight(k);
        RCv = RCv + kron(M_C(k, :)', matrixToVoigt(DW_C-lambdaCN1+wedge(lambdaGN1, CN05)+1/3*lambdacN1*GN05, 'strain')) * detJ * gaussWeight(k);
        RGv = RGv + kron(M_G(k, :)', matrixToVoigt(DW_G-lambdaGN1+1/3*lambdacN1*CN05, 'strain')) * detJ * gaussWeight(k);
        Rc = Rc + M_c(k, :)' * (DW_c - lambdacN1) * detJ * gaussWeight(k);
        RLambdaCv = RLambdaCv + kron(M_LambdaC(k, :)', matrixToVoigt(CxN1-CN1, 'strain')) * detJ * gaussWeight(k);
        RLambdaGv = RLambdaGv + kron(M_LambdaG(k, :)', matrixToVoigt(0.5*wedge(CN1, CN1)-GN1, 'strain')) * detJ * gaussWeight(k);
        RLambdac = RLambdac + M_Lambdac(k, :)' * (1 / 3 * innerProduct(GN1, CN1) - cN1) * detJ * gaussWeight(k);

        % TANGENT
        % KXX
        A1 = dN_X_I' * lambdaCN1 * dN_X_I * detJ * gaussWeight(k);
        temporaryKXX = zeros(numberOfXDofs);
        for g = 1:dimension
            temporaryKXX(g:dimension:numberOfXDofs, g:dimension:numberOfXDofs) = A1;
        end
        KXX = KXX + temporaryKXX;

        % KXLambdaC
        KXLambdaC = KXLambdaC + 2 * kron(M_LambdaC(k, :), BN05') * detJ * gaussWeight(k);

        % =====

        % KCC
        KCC = KCC + 1 / 2 * kron(M_C(k, :)'*M_C(k, :), tangentOperatorD(lambdaGN1)) * detJ * gaussWeight(k);

        % KCG
        KCG = KCG + 1 / 6 * kron(M_C(k, :)'*M_G(k, :), lambdacN1*I6) * detJ * gaussWeight(k);

        % KCLambdaC
        KCLambdaC = KCLambdaC - kron(M_C(k, :)'*M_LambdaC(k, :), I6) * detJ * gaussWeight(k);

        % KCLambdaG
        KCLambdaG = KCLambdaG + kron(M_C(k, :)'*M_LambdaG(k, :), tangentOperatorD(CN05)) * detJ * gaussWeight(k);

        % KCLambdac
        KCLambdac = KCLambdac + 1 / 3 * kron(M_C(k, :)'*M_Lambdac(k, :), matrixToVoigt(GN05, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KGC
        KGC = KGC + 1 / 6 * kron(M_G(k, :)'*M_C(k, :), lambdacN1*I6) * detJ * gaussWeight(k);

        % KGLambdaG
        KGLambdaG = KGLambdaG - kron(M_G(k, :)'*M_LambdaG(k, :), I6) * detJ * gaussWeight(k);

        % KGLambdac
        KGLambdac = KGLambdac + 1 / 3 * kron(M_G(k, :)'*M_Lambdac(k, :), matrixToVoigt(CN05, 'strain')) * detJ * gaussWeight(k);

        % =====

        % Kcc
        if flagDiscreteGradient(k, 1)
            D2W_c_cN1 = 1 / 2 * ((d / cN1 - (c * (sqrt(cN1) - 1)) / sqrt(cN1)) / (cN - cN1) + (c * (sqrt(cN) - 1)^2 - c * (sqrt(cN1) - 1)^2 - 2 * d * log(sqrt(cN)) + 2 * d * log(sqrt(cN1))) / (cN - cN1)^2);
        else
            D2W_c_cN1 = 1 / 2 * (d / (2 * cN05^2) + c / (4 * (cN05)^(3 / 2)));
        end
        Kcc = Kcc + M_c(k, :)' * M_c(k, :) * D2W_c_cN1 * detJ * gaussWeight(k);

        % KcLambdac
        KcLambdac = KcLambdac - M_c(k, :)' * M_Lambdac(k, :) * detJ * gaussWeight(k);

        % =====

        % KLambdaCX
        KLambdaCX = KLambdaCX + 2 * kron(M_LambdaC(k, :)', BN1) * detJ * gaussWeight(k);

        % KLambdaCC
        KLambdaCC = KLambdaCC - kron(M_LambdaC(k, :)'*M_C(k, :), I6) * detJ * gaussWeight(k);

        % =====

        % KLambdaGC
        KLambdaGC = KLambdaGC + kron(M_LambdaG(k, :)'*M_C(k, :), tangentOperatorD(CN1)) * detJ * gaussWeight(k);

        % KLambdaGG
        KLambdaGG = KLambdaGG - kron(M_LambdaG(k, :)'*M_G(k, :), I6) * detJ * gaussWeight(k);

        % =====

        % KLambdacC
        KLambdacC = KLambdacC + 1 / 3 * kron(M_Lambdac(k, :)'*M_C(k, :), matrixToVoigt(GN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacG
        KLambdacG = KLambdacG + 1 / 3 * kron(M_Lambdac(k, :)'*M_G(k, :), matrixToVoigt(CN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacc
        KLambdacc = KLambdacc - M_Lambdac(k, :)' * M_c(k, :) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        SN1 = 2 * lambdaCN1;
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RCv;
    rData{3} = RGv;
    rData{4} = Rc;
    rData{5} = RLambdaCv;
    rData{6} = RLambdaGv;
    rData{7} = RLambdac;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 5} = KXLambdaC;
    kData{2, 2} = KCC;
    kData{2, 3} = KCG;
    kData{2, 5} = KCLambdaC;
    kData{2, 6} = KCLambdaG;
    kData{2, 7} = KCLambdac;
    kData{3, 2} = KGC;
    kData{3, 6} = KGLambdaG;
    kData{3, 7} = KGLambdac;
    kData{4, 4} = Kcc;
    kData{4, 7} = KcLambdac;
    kData{5, 1} = KLambdaCX;
    kData{5, 2} = KLambdaCC;
    kData{6, 2} = KLambdaGC;
    kData{6, 3} = KLambdaGG;
    kData{7, 1} = KLambdacX;
    kData{7, 2} = KLambdacC;
    kData{7, 3} = KLambdacG;
    kData{7, 4} = KLambdacc;

    % pass flagDiscreteGradient
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
end
end

function D = tangentOperatorD(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -2 * A(2, 3);
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -2 * A(1, 3);
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -2 * A(1, 2);
D(4, 3) = -2 * A(1, 2);
D(4, 4) = -2 * A(3, 3);
D(4, 5) = 2 * A(1, 3);
D(4, 6) = 2 * A(2, 3);
D(5, 1) = -2 * A(2, 3);
D(5, 4) = 2 * A(1, 3);
D(5, 5) = -2 * A(1, 1);
D(5, 6) = 2 * A(1, 2);
D(6, 2) = -2 * A(1, 3);
D(6, 4) = 2 * A(2, 3);
D(6, 5) = 2 * A(1, 2);
D(6, 6) = -2 * A(2, 2);
end