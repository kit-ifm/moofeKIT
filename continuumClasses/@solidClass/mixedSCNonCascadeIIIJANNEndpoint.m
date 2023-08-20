function [rData, kData, elementEnergy, array] = mixedSCNonCascadeIIIJANNEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDSCMOONEYRIVLINENDPOINT Element routine of class solidClass.
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
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
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
% mixedSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% load trained ann model
ANNObject = obj.artificialNeuralNetworkObject;

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

Imat = eye(3);

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edAlphaN1 = dofs.edAlphaN1.';

extractedIN1 = edAlphaN1(1:numberOfInternalNodes);
extractedIIN1 = edAlphaN1(1*numberOfInternalNodes+1:2*numberOfInternalNodes);
extractedJN1 = edAlphaN1(2*numberOfInternalNodes+1:3*numberOfInternalNodes);
extractedLambdaIN1 = edAlphaN1(3*numberOfInternalNodes+1:4*numberOfInternalNodes);
extractedLambdaIIN1 = edAlphaN1(4*numberOfInternalNodes+1:5*numberOfInternalNodes);
extractedLambdaJN1 = edAlphaN1(5*numberOfInternalNodes+1:6*numberOfInternalNodes);

% initialize residual
RD = rData{1};
RI = rData{2};
RII = rData{3};
RIII = rData{4};
RLambdaI = rData{5};
RLambdaII = rData{6};
RLambdaIII = rData{7};

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
    [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    % compute the values of the variables at the current Gauss point
    IN1 = extractedIN1.' * M_k_I(k, :)';
    IIN1 = extractedIIN1.' * M_k_I(k, :)';
    JN1 = extractedJN1.' * M_k_I(k, :)';
    lambdaIN1 = extractedLambdaIN1.' * M_k_I(k, :)';
    lambdaIIN1 = extractedLambdaIIN1.' * M_k_I(k, :)';
    lambdaJN1 = extractedLambdaJN1.' * M_k_I(k, :)';

    % deformation gradient - displacement-based
    FxN1 = edN1 * dN_X_I';

    % strain measures - displacement-based
    CxN1 = FxN1.' * FxN1;

    % strain measures - displacement-based
    GxN1 = 0.5*wedge(CxN1, CxN1);

    % strain measures - displacement-based
    % cxN1 = 1/3*innerProduct(GxN1, CxN1);
    JxN1 = det(FxN1);

    % nodal operator matrix
    BN1 = BMatrix(dN_X_I, FxN1);

    % 2nd PK stress
    SN1 = 2 * (lambdaIN1*Imat + lambdaIIN1*wedge(Imat,CxN1) + 1/2*lambdaJN1*JxN1^(-1)*GxN1);

    if ~computePostData
        % ENERGY
        WN1 = ANNObject.computeEnergyANN(CxN1,GxN1,JN1^2);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + WN1 * detJ * gaussWeight(k);

        % DERIVATIVES
        [dW_I, dW_II, dW_J] = ANNObject.computeDiffEnergyANNIIIJ(IN1,IIN1,JN1);

        % RESIDUAL
        RD = RD + BN1.' * matrixToVoigt(SN1,'stress') * detJ * gaussWeight(k);
        RI = RI + M_k_I(k, :)' * (dW_I - lambdaIN1) * detJ * gaussWeight(k);
        RII = RII + M_k_I(k, :)' * (dW_II - lambdaIIN1) * detJ * gaussWeight(k);
        RIII = RIII + M_k_I(k, :)' * (dW_J - lambdaJN1) * detJ * gaussWeight(k);
        RLambdaI = RLambdaI + M_k_I(k, :)' * (trace(CxN1) - IN1) * detJ * gaussWeight(k);
        RLambdaII = RLambdaII + M_k_I(k, :)' * (trace(GxN1) - IIN1) * detJ * gaussWeight(k);
        RLambdaIII = RLambdaIII + M_k_I(k, :)' * (JxN1 - JN1) * detJ * gaussWeight(k);

        if 0
            % TANGENT
            % KXX
            A1 = 2 * dN_X_I' * lambdaCN1 * dN_X_I * detJ * gaussWeight(k);
            temporaryKXX = zeros(numberOfXDofs);
            for g = 1:dimension
                temporaryKXX(g:dimension:numberOfXDofs, g:dimension:numberOfXDofs) = A1;
            end
            KXX = KXX + temporaryKXX;

            % KXLambdaC
            KXLambdaC = KXLambdaC + 2 * kron(M_k_I(k, :), BN1') * detJ * gaussWeight(k);

            % =====

            % KCC
            KCC = KCC + kron(M_k_I(k, :)' * M_k_I(k, :), tangentOperatorD(lambdaGN1)) * detJ * gaussWeight(k);

            % KCG
            KCG = KCG + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), lambdaJN1*I6) * detJ * gaussWeight(k);

            % KCLambdaC
            KCLambdaC = KCLambdaC - kron(M_k_I(k, :)' * M_k_I(k, :), I6) * detJ * gaussWeight(k);

            % KCLambdaG
            KCLambdaG = KCLambdaG + kron(M_k_I(k, :)' * M_k_I(k, :), tangentOperatorD(CN1)) * detJ * gaussWeight(k);

            % KCLambdac
            KCLambdac = KCLambdac + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), matrixToVoigt(GN1, 'strain')) * detJ * gaussWeight(k);

            % =====

            % KGC
            KGC = KGC + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), lambdaJN1*I6) * detJ * gaussWeight(k);

            % KGLambdaG
            KGLambdaG = KGLambdaG - kron(M_k_I(k, :)' * M_k_I(k, :), I6) * detJ * gaussWeight(k);

            % KGLambdac
            KGLambdac = KGLambdac + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), matrixToVoigt(CN1, 'strain')) * detJ * gaussWeight(k);

            % =====

            % Kcc
            D2W_c_cN1 = d / (2 * IIIN1^2) + c / (4 * (IIIN1)^(3 / 2));
            Kcc = Kcc + M_k_I(k, :)' * M_k_I(k, :) * D2W_c_cN1 * detJ * gaussWeight(k);

            % KcLambdac
            KcLambdac = KcLambdac - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);

            % =====

            % KLambdaCX
            KLambdaCX = KLambdaCX + 2 * kron(M_k_I(k, :)', BN1) * detJ * gaussWeight(k);

            % KLambdaCC
            KLambdaCC = KLambdaCC - kron(M_k_I(k, :)' * M_k_I(k, :), I6) * detJ * gaussWeight(k);

            % =====

            % KLambdaGC
            KLambdaGC = KLambdaGC + kron(M_k_I(k, :)' * M_k_I(k, :), tangentOperatorD(CN1)) * detJ * gaussWeight(k);

            % KLambdaGG
            KLambdaGG = KLambdaGG - kron(M_k_I(k, :)' * M_k_I(k, :), I6) * detJ * gaussWeight(k);

            % =====

            % KLambdacC
            KLambdacC = KLambdacC + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), matrixToVoigt(GN1, 'strain')') * detJ * gaussWeight(k);

            % KLambdacG
            KLambdacG = KLambdacG + 1 / 3 * kron(M_k_I(k, :)' * M_k_I(k, :), matrixToVoigt(CN1, 'strain')') * detJ * gaussWeight(k);

            % KLambdacc
            KLambdacc = KLambdacc - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);
        end
    else
        % STRESS COMPUTATION
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RD;
    rData{2} = RI;
    rData{3} = RII;
    rData{4} = RIII;
    rData{5} = RLambdaI;
    rData{6} = RLambdaII;
    rData{7} = RLambdaIII;

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