function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDSCMOONEYRIVLINMIDPOINT Element routine of class solidClass.
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
% The routine is suitable for dynamic simulation where the midpoint rule is
% used ('Midpoint').
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
% mixedSCMooneyRivlinEndpoint
% mixedSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Felix Zaehringer, Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof;
numberOfXDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalNodes = size(M_k_I, 2);

dimension = obj.dimension;

I6 = [eye(3), zeros(3); zeros(3), 2 * eye(3)];

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edRef = obj.qR(edof(e, :), 1:dimension).';

edN = obj.qN(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
edN05 = 1 / 2 * (edN + edN1);

edAlphaN = mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedCNv = edAlphaN(1:6*numberOfInternalNodes);
extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes);
extractedGNv = edAlphaN(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedGN1v = edAlphaN1(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedcN = edAlphaN(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedcN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedLambdaCN1v = edAlphaN1(13*numberOfInternalNodes+1:19*numberOfInternalNodes);
extractedLambdaGN1v = edAlphaN1(19*numberOfInternalNodes+1:25*numberOfInternalNodes);
extractedLambdacN1 = edAlphaN1(25*numberOfInternalNodes+1:26*numberOfInternalNodes);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(edRef, dN_xi_k_I);

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

%% GAUSS LOOP
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

    cN = extractedcN.' * M_k_I(k, :)';
    cN1 = extractedcN1.' * M_k_I(k, :)';
    cN05 = 0.5 * (cN + cN1);

    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M_k_I(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');

    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M_k_I(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');

    lambdacN1 = extractedLambdacN1.' * M_k_I(k, :)';

    % deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN05 = edN05 * dN_X_I';

    % strain measures
    CxN05 = FxN05.' * FxN05;

    % nodal operator matrix
    BN05 = BMatrix(dN_X_I, FxN05);

    if ~computePostData
        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2) * detJ * gaussWeight(k);

        % DERIVATIVES
        DW_C = a * eye(3);
        DW_G = b * eye(3);
        DW_c = -d / (2 * cN05) + c / 2 * (1 - 1 / sqrt(cN05));

        % RESIDUAL
        RX = RX + 2 * BN05.' * lambdaCN1v * detJ * gaussWeight(k);
        RCv = RCv + kron(M_k_I(k, :)', matrixToVoigt(DW_C-lambdaCN1+wedge(lambdaGN1, CN05)+1/3*lambdacN1*GN05, 'stress')) * detJ * gaussWeight(k);
        RGv = RGv + kron(M_k_I(k, :)', matrixToVoigt(DW_G-lambdaGN1+1/3*lambdacN1*CN05, 'stress')) * detJ * gaussWeight(k);
        Rc = Rc + M_k_I(k, :)' * (DW_c - lambdacN1) * detJ * gaussWeight(k);
        RLambdaCv = RLambdaCv + kron(M_k_I(k, :)', matrixToVoigt(CxN05-CN05, 'stress')) * detJ * gaussWeight(k);
        RLambdaGv = RLambdaGv + kron(M_k_I(k, :)', matrixToVoigt(0.5*wedge(CN05, CN05)-GN05, 'stress')) * detJ * gaussWeight(k);
        RLambdac = RLambdac + M_k_I(k, :)' * (1 / 3 * innerProduct(GN05, CN05) - cN05) * detJ * gaussWeight(k);

        % TANGENT
        % KXX
        A1 = dN_X_I' * lambdaCN1 * dN_X_I * detJ * gaussWeight(k);
        temporaryKXX = zeros(numberOfXDofs);
        for g = 1:dimension
            temporaryKXX(g:dimension:numberOfXDofs, g:dimension:numberOfXDofs) = A1;
        end
        KXX = KXX + temporaryKXX;

        % KXLambdaC
        KXLambdaC = KXLambdaC + 2 * kron(M_k_I(k, :), BN05') * detJ * gaussWeight(k);

        % =====

        % KCC
        KCC = KCC + 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(lambdaGN1)) * detJ * gaussWeight(k);

        % KCG
        KCG = KCG + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*I6) * detJ * gaussWeight(k);

        % KCLambdaC
        KCLambdaC = KCLambdaC - kron(M_k_I(k, :)'*M_k_I(k, :), I6) * detJ * gaussWeight(k);

        % KCLambdaG
        KCLambdaG = KCLambdaG + kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(CN05)) * detJ * gaussWeight(k);

        % KCLambdac
        KCLambdac = KCLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN05, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KGC
        KGC = KGC + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*I6) * detJ * gaussWeight(k);

        % KGLambdaG
        KGLambdaG = KGLambdaG - kron(M_k_I(k, :)'*M_k_I(k, :), I6) * detJ * gaussWeight(k);

        % KGLambdac
        KGLambdac = KGLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN05, 'strain')) * detJ * gaussWeight(k);

        % =====

        % Kcc
        D2W_c_cN1 = 1 / 2 * (d / (2 * cN05^2) + c / (4 * (cN05)^(3 / 2)));
        Kcc = Kcc + M_k_I(k, :)' * M_k_I(k, :) * D2W_c_cN1 * detJ * gaussWeight(k);

        % KcLambdac
        KcLambdac = KcLambdac - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);

        % =====

        % KLambdaCX
        KLambdaCX = KLambdaCX + kron(M_k_I(k, :)', BN05) * detJ * gaussWeight(k);

        % KLambdaCC
        KLambdaCC = KLambdaCC - 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), I6) * detJ * gaussWeight(k);

        % =====

        % KLambdaGC
        KLambdaGC = KLambdaGC + 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(CN05)) * detJ * gaussWeight(k);

        % KLambdaGG
        KLambdaGG = KLambdaGG - 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), I6) * detJ * gaussWeight(k);

        % =====

        % KLambdacC
        KLambdacC = KLambdacC + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN05, 'strain')') * detJ * gaussWeight(k);

        % KLambdacG
        KLambdacG = KLambdacG + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN05, 'strain')') * detJ * gaussWeight(k);

        % KLambdacc
        KLambdacc = KLambdacc - 1 / 2 * M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        SN1 = 2 * lambdaCN1;
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
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