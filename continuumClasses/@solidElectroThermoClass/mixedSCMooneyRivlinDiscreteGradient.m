function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDSCMOONEYRIVLINDISCRETEGRADIENT Element routine of class solidElectroThermoClass.
%
% FORMULATION
% This is a 'mixed'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic  Mooney-Rivlin
% ('MooneyRivlin') model (nonlinear geometric and nonlinear material/
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
% mixedSCMooneyRivlinDiscreteGradient(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, P, T, D, C, G, c, LambdaC, LambdaG, Lambdac)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, P, T, D, C, G, c, LambdaC,
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
% Master thesis Felix Zaehringer, 2021
%
% SEE ALSO
% mixedSCMooneyRivlinMidpoint
% mixedSCMooneyRivlinEndpoint
%
% CREATOR(S)
% Felix Zaehringer, Marlon Franke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;
% element degree of freedom tables and more
edof = meshObject.edof(e, :);
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
% material data and voigt notation
numberOfInternalNodes = size(M_k_I, 2);
dimension = obj.dimension;
DT = setupObject.timeStepSize;

% map
map.Voigt = [1, 5, 9, 4, 8, 7]';
map.VoigtInv = [1, 4, 6; 4, 2, 5; 6, 5, 3];
map.VoigtFull = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
map.VoigtFullInv = [1, 4, 6; 7, 2, 5; 9, 8, 3];
map.Isym = [eye(3), zeros(3); zeros(3), 2 * eye(3)];
map.Isyminv = [eye(3), zeros(3); zeros(3), 1 / 2 * eye(3)];

% initialize energy
elementEnergy.internalEnergy = 0;
elementEnergy.helmholtz = 0;
elementEnergy.DE = 0;
elementEnergy.TS = 0;
elementEnergy.S = 0;
elementEnergy.deltaS = 0;
elementEnergy.deltaU = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 5)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

%% extract dofs for element
edN1 = dofs.edN1;
phiN1 = dofs.phiN1;
thetaN1 = dofs.thetaN1;
edAlphaN1 = dofs.edAlphaN1;

% displacement DOFs
edN = obj.qN(edof, 1:dimension)';
edN05 = 1 / 2 * (edN + edN1);
phiN = obj.qN(edof, dimension+1)';
phiN05 = 1 / 2 * (phiN + phiN1);
edAlphaN = obj.mixedFEObject.qN(e, :);
thetaN = obj.qN(edof, dimension+2)';
thetaN05 = 1 / 2 * (thetaN1 + thetaN);

% internal DOFs
extractedDN = edAlphaN(1:3*numberOfInternalNodes).';
extractedCNv = edAlphaN(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
extractedGNv = edAlphaN(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
extractedcN = edAlphaN(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';

extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';
extractedCN1v = edAlphaN1(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
extractedGN1v = edAlphaN1(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
extractedcN1 = edAlphaN1(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
extractedLambdaCN1v = edAlphaN1(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
extractedLambdaGN1v = edAlphaN1(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
extractedLambdacN1 = edAlphaN1(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';

%
edR = obj.qR(edof, 1:dimension)';

%% Material parameters
a = materialObject.a;
b = materialObject.b;
c1 = materialObject.c1;
d1 = materialObject.d1;
c2 = materialObject.c2;
% d2 = materialObject.d2;
e0 = materialObject.e0;
er = materialObject.e1;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;
% rho0 = materialObject.rhoSource;
% R = materialObject.RSource;

%% Initialization
RX = rData{1, 1};
RP = rData{2, 1};
RT = rData{3, 1};
RD = rData{4, 1};
RCv = rData{5, 1};
RGv = rData{6, 1};
Rc = rData{7, 1};
RLambdaCv = rData{8, 1};
RLambdaGv = rData{9, 1};
RLambdac = rData{10, 1};
KXX=kData{1,1}; KXP=kData{1,2}; KXT=kData{1,3}; KXD=kData{1,4}; KXC=kData{1,5}; KXG=kData{1,6}; KXc=kData{1,7}; KXLambdaC=kData{1,8}; KXLambdaG=kData{1,9}; KXLambdac=kData{1,10}; 
KPX=kData{2,1}; KPP=kData{2,2}; KPT=kData{2,3}; KPD=kData{2,4}; KPC=kData{2,5}; KPG=kData{2,6}; KPc=kData{2,7}; KPLambdaC=kData{2,8}; KPLambdaG=kData{2,9}; KPLambdac=kData{2,10}; 
KTX=kData{3,1}; KTP=kData{3,2}; KTT=kData{3,3}; KTD=kData{3,4}; KTC=kData{3,5}; KTG=kData{3,6}; KTc=kData{3,7}; KTLambdaC=kData{3,8}; KTLambdaG=kData{3,9}; KTLambdac=kData{3,10}; 
KDX=kData{4,1}; KDP=kData{4,2}; KDT=kData{4,3}; KDD=kData{4,4}; KDC=kData{4,5}; KDG=kData{4,6}; KDc=kData{4,7}; KDLambdaC=kData{4,8}; KDLambdaG=kData{4,9}; KDLambdac=kData{4,10}; 
KCX=kData{5,1}; KCP=kData{5,2}; KCT=kData{5,3}; KCD=kData{5,4}; KCC=kData{5,5}; KCG=kData{5,6}; KCc=kData{5,7}; KCLambdaC=kData{5,8}; KCLambdaG=kData{5,9}; KCLambdac=kData{5,10}; 
KGX=kData{6,1}; KGP=kData{6,2}; KGT=kData{6,3}; KGD=kData{6,4}; KGC=kData{6,5}; KGG=kData{6,6}; KGc=kData{6,7}; KGLambdaC=kData{6,8}; KGLambdaG=kData{6,9}; KGLambdac=kData{6,10}; 
KcX=kData{7,1}; KcP=kData{7,2}; KcT=kData{7,3}; KcD=kData{7,4}; KcC=kData{7,5}; KcG=kData{7,6}; Kcc=kData{7,7}; KcLambdaC=kData{7,8}; KcLambdaG=kData{7,9}; KcLambdac=kData{7,10}; 
KLambdaCX=kData{8,1}; KLambdaCP=kData{8,2}; KLambdaCT=kData{8,3}; KLambdaCD=kData{8,4}; KLambdaCC=kData{8,5}; KLambdaCG=kData{8,6}; KLambdaCc=kData{8,7}; KLambdaCLambdaC=kData{8,8}; KLambdaCLambdaG=kData{8,9}; KLambdaCLambdac=kData{8,10}; 
KLambdaGX=kData{9,1}; KLambdaGP=kData{9,2}; KLambdaGT=kData{9,3}; KLambdaGD=kData{9,4}; KLambdaGC=kData{9,5}; KLambdaGG=kData{9,6}; KLambdaGc=kData{9,7}; KLambdaGLambdaC=kData{9,8}; KLambdaGLambdaG=kData{9,9}; KLambdaGLambdac=kData{9,10}; 
KLambdacX=kData{10,1}; KLambdacP=kData{10,2}; KLambdacT=kData{10,3}; KLambdacD=kData{10,4}; KLambdacC=kData{10,5}; KLambdacG=kData{10,6}; KLambdacc=kData{10,7}; KLambdacLambdaC=kData{10,8}; KLambdacLambdaG=kData{10,9}; KLambdacLambdac=kData{10,10};

%% DOFs
numberOfXDOFs = numel(RX);

%% Gauss loop
for k = 1:numberOfGausspoints
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);

    % assembly internal dofs
    DN = reshape(extractedDN, 3, []) * M_k_I(k, :)';
    CNv = reshape(extractedCNv, 6, []) * M_k_I(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    GNv = reshape(extractedGNv, 6, []) * M_k_I(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    cN = extractedcN' * M_k_I(k, :)';

    DN1 = reshape(extractedDN1, 3, []) * M_k_I(k, :)';
    CN1v = reshape(extractedCN1v, 6, []) * M_k_I(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M_k_I(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    cN1 = extractedcN1.' * M_k_I(k, :)';
    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M_k_I(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');
    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M_k_I(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');
    lambdacN1 = extractedLambdacN1.' * M_k_I(k, :)';

    DN05 = 0.5 * (DN + DN1);
    CN05 = 0.5 * (CN + CN1);
    GN05 = 0.5 * (GN + GN1);
    cN05 = 0.5 * (cN + cN1);

    % temperature-based quantities
    thetaN1e = N_k_I(k, :) * thetaN1.';
    thetaNe = N_k_I(k, :) * thetaN.';
    thetaN05e = 0.5 * (thetaN1e+thetaNe);
    etaN = dimension * beta * c2 * (cN - 1) + kappa * log(thetaNe/thetaR);
    etaN1 = dimension * beta * c2 * (cN1 - 1) + kappa * log(thetaN1e/thetaR);

    % deformation-based quantities
    FxN1 = edN1 * dN_X_I';
    FxN = edN * dN_X_I';
    FxN05 = edN05 * dN_X_I';
    CxN1 = FxN1.' * FxN1;
    CxN = FxN.' * FxN;
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    GxN = 0.5 * wedge(CxN, CxN);
    cxN1 = det(CxN1);
    cxN = det(CxN);

    % electrical quantities
    EN1 = -dN_X_I * phiN1.';
    EN = -dN_X_I * phiN.';

    % B-matrix (midpoint configuration)
    BN05 = BMatrix(dN_X_I, FxN05);

    % B-matrix
    BN1 = BMatrix(dN_X_I, FxN1);

    %% Stress computation
    if computePostData
        SN1 = 2 * lambdaCN1;
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        stressTensor.D = DN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end

    %% energy
    PsiTN = kappa * (thetaNe - thetaR - thetaNe * log(thetaNe/thetaR));
    PsiTMN = -dimension * beta * c2 * (thetaNe - thetaR) * (cN - 1);
    PsiEMN = 1 / (2 * er * e0 * cN^(1 / 2)) * DN' * CN * DN;
    PsiMIsoN = a * (trace(CN) - 3) + b * (trace(GN) - 3);
    PsiMVolN = c1 / 2 * (sqrt(cN) - 1)^2 - d1 * log(sqrt(cN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
    u0N = PsiN - DN' * EN + thetaNe * etaN;
    PsiTN1 = kappa * (thetaN1e-thetaR - thetaN1e * log(thetaN1e/thetaR));
    PsiTMN1 = -dimension * beta * c2 * (thetaN1e-thetaR) * (cN1 - 1);
    PsiEMN1 = 1 / (2 * er * e0 * cN1^(1 / 2)) * DN1' * CN1 * DN1;
    PsiMIsoN1 = a * (trace(CN1) - 3) + b * (trace(GN1) - 3);
    PsiMVolN1 = c1 / 2 * (sqrt(cN1) - 1)^2 - d1 * log(sqrt(cN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
    u0N1 = PsiN1 - DN1' * EN1 + thetaN1e * etaN1;

    elementEnergy.internalEnergy = elementEnergy.internalEnergy + u0N1 * detJ * gaussWeight(k);
    elementEnergy.helmholtz = elementEnergy.helmholtz + PsiN1 * detJ * gaussWeight(k);
    elementEnergy.DE = elementEnergy.DE + DN1' * EN1 * detJ * gaussWeight(k);
    elementEnergy.TS = elementEnergy.TS + thetaN1e * etaN1 * detJ * gaussWeight(k);
    elementEnergy.S = elementEnergy.S + etaN1 * detJ * gaussWeight(k);
    elementEnergy.deltaS = elementEnergy.deltaS + (etaN1 - etaN) * detJ * gaussWeight(k);
    elementEnergy.deltaU = elementEnergy.deltaU + (u0N1 - u0N) * detJ * gaussWeight(k);

    %% Discrete Gradient
    deltaC = CN1 - CN;
    normDeltaC = sqrt(innerProduct(deltaC, deltaC));
    deltac = cN1 - cN;
    normDeltac = abs(deltac);
    deltaD = DN1 - DN;
    normDeltaD = sqrt(deltaD'*deltaD);
    deltaT = thetaN1e-thetaNe;
    normDeltaT = abs(deltaT);

    % conditions
    if ~flagNumericalTangent
        discreteGradientCondition_T = (normDeltaT > 10 * eps);
        discreteGradientCondition_D = (normDeltaD > 10 * eps);
        discreteGradientCondition_C = (normDeltaC > 10 * eps);
        discreteGradientCondition_c = (normDeltac > 10 * eps);
    else
        discreteGradientCondition_T = flagDiscreteGradient(k, 1);
        discreteGradientCondition_D = flagDiscreteGradient(k, 2);
        discreteGradientCondition_C = flagDiscreteGradient(k, 3);
        discreteGradientCondition_c = flagDiscreteGradient(k, 4);
    end

    % discrete derivatives
    if discreteGradientCondition_T
        flagDiscreteGradient(k, 1) = 1;
        DW_theta = -3 / 2 * beta * c2 * (cN1 + cN - 2) ...
            +kappa / deltaT * (thetaN1e * (1 - log(thetaN1e/thetaR)) - thetaNe * (1 - log(thetaNe/thetaR)));
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_theta = -dimension * beta * c2 * (cN05 - 1) - kappa * log(thetaN05e/thetaR);
    end

    if discreteGradientCondition_D
        flagDiscreteGradient(k, 2) = 1;
        DW_D = 1 / (2 * er * e0) * (CN1 / sqrt(cN1) + CN / sqrt(cN)) * DN05;
    else
        flagDiscreteGradient(k, 2) = 0;
        DW_D = 1 / (er * e0 * sqrt(cN05)) * CN05 * DN05;
    end

    if discreteGradientCondition_C
        flagDiscreteGradient(k, 3) = 1;
        DW_C = a * eye(3) + 1 / (4 * er * e0) * (1 / sqrt(cN1) * (DN1 * DN1.') + 1 / sqrt(cN) * (DN * DN.'));
    else
        flagDiscreteGradient(k, 3) = 0;
        DW_C = a * eye(3) + 1 / (2 * er * e0 * sqrt(cN05)) * (DN05 * DN05.');
    end

    DW_G = b * eye(3);

    if discreteGradientCondition_c
        flagDiscreteGradient(k, 4) = 1;
        DW_c = d1 / deltac * log(sqrt(cN/cN1)) ...
            +c1 / (2 * deltac) * ((sqrt(cN1) - 1)^2 - (sqrt(cN) - 1)^2) ...
            +1 / (4 * er * e0 * deltac) * (1 / sqrt(cN1) - 1 / sqrt(cN)) * (DN.' * (CN1 * DN) + DN1.' * (CN * DN1)) ...
            -3 / 2 * beta * c2 * (thetaN1e+thetaNe - 2 * thetaR);
    else
        flagDiscreteGradient(k, 4) = 0;
        DW_c = -d1 / (2 * cN05) + c1 / 2 * (1 - 1 / sqrt(cN05)) ...
            -1 / (4 * er * e0 * cN05^(3 / 2)) * DN05.' * (CN05 * DN05) ...
            -dimension * beta * c2 * (thetaN05e-thetaR);
    end

    %% Residual
    etaAlgo = -DW_theta;
    QN05 = -k0 * cN05^(-1) * GN05 * (dN_X_I * thetaN05.');

    RX = RX + BN05.' * 2 * matrixToVoigt(lambdaCN1, 'stress') * detJ * gaussWeight(k);
    RT = RT + (N_k_I(k, :)' * (thetaN1e * etaN1 - thetaNe * etaN) / DT + N_k_I(k, :)' * (thetaN1e-thetaNe) / DT * DW_theta - dN_X_I' * QN05) * detJ * gaussWeight(k);
    RP = RP + (dN_X_I' * DN05) * detJ * gaussWeight(k);
    RD = RD + kron(M_k_I(k, :)', (DW_D + dN_X_I * phiN05.')) * detJ * gaussWeight(k);
    RCv = RCv + kron(M_k_I(k, :)', matrixToVoigt(DW_C-lambdaCN1+wedge(lambdaGN1, CN05)+1/3*lambdacN1*GN05, 'strain')) * detJ * gaussWeight(k);
    RGv = RGv + kron(M_k_I(k, :)', matrixToVoigt(DW_G-lambdaGN1+1/3*lambdacN1*CN05, 'strain')) * detJ * gaussWeight(k);
    Rc = Rc + M_k_I(k, :)' * (DW_c - lambdacN1) * detJ * gaussWeight(k);
    RLambdaCv = RLambdaCv + kron(M_k_I(k, :)', matrixToVoigt(FxN1.'*FxN1-CN1, 'strain')) * detJ * gaussWeight(k);
    RLambdaGv = RLambdaGv + kron(M_k_I(k, :)', matrixToVoigt(0.5*wedge(CN1, CN1)-GN1, 'strain')) * detJ * gaussWeight(k);
    RLambdac = RLambdac + M_k_I(k, :)' * (innerProduct(1/3*GN1, CN1) - cN1) * detJ * gaussWeight(k);

    %% Tangent
    % KXX
    geometricalTangentPart = dN_X_I' * lambdaCN1 * dN_X_I;
    geometricalTangent = zeros(numberOfXDOFs);
    for g = 1:dimension
        geometricalTangent(g:dimension:numberOfXDOFs, g:dimension:numberOfXDOFs) = geometricalTangentPart;
    end
    KXX = KXX + geometricalTangent * detJ * gaussWeight(k);

    % KXLambdaC
    KXLambdaC = KXLambdaC + 2 * kron(M_k_I(k, :), BN05') * detJ * gaussWeight(k);

    % =====

    % KTT
    kTT1 = N_k_I(k, :)' * N_k_I(k, :) * 1 / DT * (etaN1 - etaAlgo);
    if flagDiscreteGradient(k, 1)
        kTT2 = N_k_I(k, :)' * N_k_I(k, :) * 1 / DT * kappa * (1 - log(thetaN1e/thetaR)) ...
            -N_k_I(k, :)' * N_k_I(k, :) * 1 / DT * kappa * thetaN1e / (thetaN1e-thetaNe) * (1 - log(thetaN1e/thetaR)) ...
            +N_k_I(k, :)' * N_k_I(k, :) * 1 / DT * kappa * thetaNe / (thetaN1e-thetaNe) * (1 - log(thetaNe/thetaR));
    else
        kTT2 = N_k_I(k, :)' * N_k_I(k, :) * 1 / DT * kappa * (1 - (thetaN1e-thetaNe) / (2 * thetaN05e));
    end
    kTT3 = 1 / 2 * dN_X_I' * k0 / cN05 * GN05 * dN_X_I;
    KTT = KTT + (kTT1 + kTT2 + kTT3) * detJ * gaussWeight(k);

    % KTG
    KTG = KTG + 1 / 2 * kron(M_k_I(k, :), dN_X_I'*k0/cN05*BErs(dN_X_I*thetaN05')) * detJ * gaussWeight(k);

    % KTc
    kTc1 = 3 / 2 * kron(N_k_I(k, :)'*M_k_I(k, :), 1/DT*(thetaN1e+thetaNe)*beta*c2) ...
        -1 / 2 * kron(M_k_I(k, :), dN_X_I'*k0/(cN05^2)*GN05*(dN_X_I * thetaN05'));
    KTc = KTc + kTc1 * detJ * gaussWeight(k);

    % =====

    % KPD
    KPD = KPD + 1 / 2 * kron(M_k_I(k, :), dN_X_I') * detJ * gaussWeight(k);

    % =====

    % KDP
    KDP = KDP + 1 / 2 * kron(M_k_I(k, :)', dN_X_I) * detJ * gaussWeight(k);

    % KDD
    if flagDiscreteGradient(k, 2)
        kDD1 = 1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0)*(CN1 / sqrt(cN1) + CN / sqrt(cN)));
    else
        kDD1 = 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * sqrt(cN05))*CN05);
    end
    KDD = KDD + kDD1 * detJ * gaussWeight(k);

    % KDC
    if flagDiscreteGradient(k, 2)
        kDC1 = 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * sqrt(cN1))*BErs(DN05));
    else
        kDC1 = 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * sqrt(cN05))*BErs(DN05));
    end
    KDC = KDC + kDC1 * detJ * gaussWeight(k);

    % KDc
    if flagDiscreteGradient(k, 2)
        KDc1 = -1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * (cN1)^(3 / 2))*CN1*DN05);
    else
        KDc1 = -1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * (cN05)^(3 / 2))*CN05*DN05);
    end
    KDc = KDc + KDc1 * detJ * gaussWeight(k);

    % =====

    % KCD
    if flagDiscreteGradient(k, 3)
        kCD1 = 1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * sqrt(cN1))*AErs(DN1));
    else
        kCD1 = 1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * sqrt(cN05))*AErs(DN05));
    end
    KCD = KCD + kCD1 * detJ * gaussWeight(k);

    % KCC
    KCC = KCC + 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(lambdaGN1)) * detJ * gaussWeight(k);

    % KCG
    KCG = KCG + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

    % KCc
    if flagDiscreteGradient(k, 3)
        KCc1 = -1 / 8 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * (cN1)^(3 / 2))*matrixToVoigt(DN1*DN1', 'strain'));
    else
        KCc1 = -1 / 8 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * (cN05)^(3 / 2))*matrixToVoigt(DN05*DN05', 'strain'));
    end
    KCc = KCc + KCc1 * detJ * gaussWeight(k);

    % KCLambdaC
    KCLambdaC = KCLambdaC - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

    % KCLambdaG
    KCLambdaG = KCLambdaG + kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(CN05)) * detJ * gaussWeight(k);

    % KCLambdac
    KCLambdac = KCLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN05, 'strain')) * detJ * gaussWeight(k);

    % =====

    % KGC
    KGC = KGC + 1 / 6 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

    % KGLambdaG
    KGLambdaG = KGLambdaG - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

    % KGLambdac
    KGLambdac = KGLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN05, 'strain')) * detJ * gaussWeight(k);

    % =====

    % KcT
    KcT = KcT - 3 / 2 * M_k_I(k, :)' * N_k_I(k, :) * beta * c2 * detJ * gaussWeight(k);

    % KcD
    if flagDiscreteGradient(k, 4)
        KcD1 = 1 / 2 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * deltac)*(1 / sqrt(cN1) - 1 / sqrt(cN))*(CN * DN1)');
    else
        KcD1 = -1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * cN05^(3 / 2))*(CN05 * DN05)');
    end
    KcD = KcD + KcD1 * detJ * gaussWeight(k);

    % KcC
    if flagDiscreteGradient(k, 4)
        KcC1 = 1 / 4 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * deltac)*(1 / sqrt(cN1) - 1 / sqrt(cN))*DN'*BErs(DN));
    else
        KcC1 = -1 / 8 * kron(M_k_I(k, :)'*M_k_I(k, :), 1/(er * e0 * cN05^(3 / 2))*DN05'*BErs(DN05));
    end
    KcC = KcC + KcC1 * detJ * gaussWeight(k);

    % Kcc
    if flagDiscreteGradient(k, 4)
        D2W_c_cN1 = -1 / 2 * d1 / (cN1 * deltac) ...
            -d1 / (deltac)^2 * log(sqrt(cN/cN1)) ...
            +1 / 2 * c1 / deltac * (1 - 1 / sqrt(cN1)) ...
            -1 / 2 * c1 / (deltac)^2 * ((sqrt(cN1) - 1)^2 - (sqrt(cN) - 1)^2) ...
            -1 / 8 * 1 / (er * e0 * (cN1)^(3 / 2) * deltac) * (DN' * (CN1 * DN) + DN1' * (CN * DN1)) ...
            -1 / 4 * 1 / (er * e0 * deltac^2) * (1 / sqrt(cN1) - 1 / sqrt(cN)) * (DN' * (CN1 * DN) + DN1' * (CN * DN1));

        Kcc1 = M_k_I(k, :)' * M_k_I(k, :) * D2W_c_cN1;
    else
        Kcc1 = M_k_I(k, :)' * M_k_I(k, :) * (d1 / (4 * cN05^2) + c1 / (8 * cN05^(3 / 2))) ...
            +M_k_I(k, :)' * M_k_I(k, :) * (3 / (16 * er * e0 * cN05^(5 / 2)) * DN05' * (CN05 * DN05));
    end
    Kcc = Kcc + Kcc1 * detJ * gaussWeight(k);

    % KcLambdac
    KcLambdac = KcLambdac - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);

    % =====

    % KLambdaCX
    KLambdaCX = KLambdaCX + 2 * kron(M_k_I(k, :)', BN1) * detJ * gaussWeight(k);

    % KLambdaCC
    KLambdaCC = KLambdaCC - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

    % =====

    % KLambdaGC
    KLambdaGC = KLambdaGC + kron(M_k_I(k, :)'*M_k_I(k, :), tangentOperatorD(CN1)) * detJ * gaussWeight(k);

    % KLambdaGG
    KLambdaGG = KLambdaGG - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

    % =====

    % KLambdacC
    KLambdacC = KLambdacC + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN1, 'strain')') * detJ * gaussWeight(k);

    % KLambdacG
    KLambdacG = KLambdacG + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN1, 'strain')') * detJ * gaussWeight(k);

    % KLambdacc
    KLambdacc = KLambdacc - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);
end
if ~computePostData
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
    rData{1} = RX;
    rData{2} = RP;
    rData{3} = RT;
    rData{4} = RD;
    rData{5} = RCv;
    rData{6} = RGv;
    rData{7} = Rc;
    rData{8} = RLambdaCv;
    rData{9} = RLambdaGv;
    rData{10} = RLambdac;
    kData{1, 1} = KXX;
    kData{1, 2} = KXP;
    kData{1, 3} = KXT;
    kData{1, 4} = KXD;
    kData{1, 5} = KXC;
    kData{1, 6} = KXG;
    kData{1, 7} = KXc;
    kData{1, 8} = KXLambdaC;
    kData{1, 9} = KXLambdaG;
    kData{1, 10} = KXLambdac;
    kData{2, 1} = KPX;
    kData{2, 2} = KPP;
    kData{2, 3} = KPT;
    kData{2, 4} = KPD;
    kData{2, 5} = KPC;
    kData{2, 6} = KPG;
    kData{2, 7} = KPc;
    kData{2, 8} = KPLambdaC;
    kData{2, 9} = KPLambdaG;
    kData{2, 10} = KPLambdac;
    kData{3, 1} = KTX;
    kData{3, 2} = KTP;
    kData{3, 3} = KTT;
    kData{3, 4} = KTD;
    kData{3, 5} = KTC;
    kData{3, 6} = KTG;
    kData{3, 7} = KTc;
    kData{3, 8} = KTLambdaC;
    kData{3, 9} = KTLambdaG;
    kData{3, 10} = KTLambdac;
    kData{4, 1} = KDX;
    kData{4, 2} = KDP;
    kData{4, 3} = KDT;
    kData{4, 4} = KDD;
    kData{4, 5} = KDC;
    kData{4, 6} = KDG;
    kData{4, 7} = KDc;
    kData{4, 8} = KDLambdaC;
    kData{4, 9} = KDLambdaG;
    kData{4, 10} = KDLambdac;
    kData{5, 1} = KCX;
    kData{5, 2} = KCP;
    kData{5, 3} = KCT;
    kData{5, 4} = KCD;
    kData{5, 5} = KCC;
    kData{5, 6} = KCG;
    kData{5, 7} = KCc;
    kData{5, 8} = KCLambdaC;
    kData{5, 9} = KCLambdaG;
    kData{5, 10} = KCLambdac;
    kData{6, 1} = KGX;
    kData{6, 2} = KGP;
    kData{6, 3} = KGT;
    kData{6, 4} = KGD;
    kData{6, 5} = KGC;
    kData{6, 6} = KGG;
    kData{6, 7} = KGc;
    kData{6, 8} = KGLambdaC;
    kData{6, 9} = KGLambdaG;
    kData{6, 10} = KGLambdac;
    kData{7, 1} = KcX;
    kData{7, 2} = KcP;
    kData{7, 3} = KcT;
    kData{7, 4} = KcD;
    kData{7, 5} = KcC;
    kData{7, 6} = KcG;
    kData{7, 7} = Kcc;
    kData{7, 8} = KcLambdaC;
    kData{7, 9} = KcLambdaG;
    kData{7, 10} = KcLambdac;
    kData{8, 1} = KLambdaCX;
    kData{8, 2} = KLambdaCP;
    kData{8, 3} = KLambdaCT;
    kData{8, 4} = KLambdaCD;
    kData{8, 5} = KLambdaCC;
    kData{8, 6} = KLambdaCG;
    kData{8, 7} = KLambdaCc;
    kData{8, 8} = KLambdaCLambdaC;
    kData{8, 9} = KLambdaCLambdaG;
    kData{8, 10} = KLambdaCLambdac;
    kData{9, 1} = KLambdaGX;
    kData{9, 2} = KLambdaGP;
    kData{9, 3} = KLambdaGT;
    kData{9, 4} = KLambdaGD;
    kData{9, 5} = KLambdaGC;
    kData{9, 6} = KLambdaGG;
    kData{9, 7} = KLambdaGc;
    kData{9, 8} = KLambdaGLambdaC;
    kData{9, 9} = KLambdaGLambdaG;
    kData{9, 10} = KLambdaGLambdac;
    kData{10, 1} = KLambdacX;
    kData{10, 2} = KLambdacP;
    kData{10, 3} = KLambdacT;
    kData{10, 4} = KLambdacD;
    kData{10, 5} = KLambdacC;
    kData{10, 6} = KLambdacG;
    kData{10, 7} = KLambdacc;
    kData{10, 8} = KLambdacLambdaC;
    kData{10, 9} = KLambdacLambdaG;
    kData{10, 10} = KLambdacLambdac;
end
end

function out = AErs(a)
out = zeros(6, 3);
if isvector(a)
    out = 2 * [a(1), 0, 0; 0, a(2), 0; 0, 0, a(3); a(2), a(1), 0; 0, a(3), a(2); a(3), 0, a(1)];
end
end

function out = BErs(b)
out = zeros(3, 6);
if isvector(b)
    out = [b(1), 0, 0, b(2), 0, b(3); 0, b(2), 0, b(1), b(3), 0; 0, 0, b(3), 0, b(2), b(1)];
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