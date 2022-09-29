function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% MIXEDSCMOONEYRIVLINENDPOINT Element routine of class solidClass.
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
% mixedSCMooneyRivlinMidpoint,
% mixedSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Felix Zaehringer, Marlon Franke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof(e, :);
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
M = mixedFEObject.shapeFunctionObject.N;
% material data and voigt notation
numberOfInternalNodes = size(M, 2);
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

%% extract dofs for element
edN1 = dofs.edN1;
phiN1 = dofs.phiN1;
thetaN1 = dofs.thetaN1;
edAlphaN1 = dofs.edAlphaN1;

% displacement DOFs
edN = obj.qN(edof, 1:dimension)';
phiN = obj.qN(edof, dimension+1)';
edAlphaN = obj.mixedFEObject.qN(e, :);
thetaN = obj.qN(edof, dimension+2)';

% internal DOFs
extractedDN = edAlphaN(1:3*numberOfInternalNodes).';
extractedcN = edAlphaN(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';
extractedCN1v = edAlphaN1(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
extractedGN1v = edAlphaN1(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
extractedcN1 = edAlphaN1(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
extractedLambdaCN1v = edAlphaN1(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
extractedLambdaGN1v = edAlphaN1(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
extractedLambdacN1 = edAlphaN1(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';

%
edRef = obj.qR(edof, 1:dimension)';

%% Jacobi-Matrix
J = edRef * dNr';
JN1 = edN1 * dNr';

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
    index = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);
    % assembly internal dofs
    DN = reshape(extractedDN, 3, []) * M(k, :)';
    cN = extractedcN.' * M(k, :)';
    DN1 = reshape(extractedDN1, 3, []) * M(k, :)';
    CN1v = reshape(extractedCN1v, 6, []) * M(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    cN1 = extractedcN1.' * M(k, :)';
    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');
    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');
    lambdacN1 = extractedLambdacN1.' * M(k, :)';
    % Temperature
    thetaN1e = N(k, :) * thetaN1.';
    thetaNe = N(k, :) * thetaN.';
    etaN = dimension * beta * c2 * (cN - 1) + kappa * log(thetaNe/thetaR);
    etaN1 = dimension * beta * c2 * (cN1 - 1) + kappa * log(thetaN1e/thetaR);
    % Deformation Gradient
    FxN1 = edN1 * dNx';
    CxN1 = FxN1.' * FxN1;
    % Electrical quantities
    EN1 = -dNx * phiN1.';
    EN = -dNx * phiN.';
    % B-matrix (midpoint configuration)
    BN1 = BMatrix(dNx, FxN1);

    if ~computePostData

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

        %% Derivatives
        DW_D = 1 / (er * e0 * sqrt(cN1)) * CN1 * DN1;
        DW_theta = -dimension * beta * c2 * (cN1 - 1) - kappa * log(thetaN1e/thetaR);
        DW_C = a * eye(3) + 1 / (2 * er * e0 * sqrt(cN1)) * (DN1 * DN1.');
        DW_G = b * eye(3);
        DW_c = -d1 / (2 * cN1) + c1 / 2 * (1 - 1 / sqrt(cN1)) - 1 / (4 * er * e0 * cN1^(3 / 2)) * DN1.' * (CN1 * DN1) - dimension * beta * c2 * (thetaN1e-thetaR);

        %% Constitutive equations
        etaN1 = dimension * beta * c2 * (cN1 - 1) + kappa * log(thetaN1e/thetaR);
        QN1 = -k0 * cN1^(-1) * GN1 * (dNx * thetaN1.');

        %% Residual
        RX = RX + BN1.' * 2 * matrixToVoigt(lambdaCN1, 'stress') * detJ * gaussWeight(k);
        RT = RT + (N(k, :)' * (thetaN1e * etaN1 - thetaNe * etaN) / DT + N(k, :)' * (thetaN1e-thetaNe) / DT * DW_theta - dNx' * QN1) * detJ * gaussWeight(k);
        RP = RP + (dNx' * DN1) * detJ * gaussWeight(k);
        RD = RD + kron(M(k, :)', (DW_D + dNx * phiN1.')) * detJ * gaussWeight(k);
        RCv = RCv + kron(M(k, :)', matrixToVoigt(DW_C-lambdaCN1+wedge(lambdaGN1, CN1)+1/3*lambdacN1*GN1, 'strain')) * detJ * gaussWeight(k);
        RGv = RGv + kron(M(k, :)', matrixToVoigt(DW_G-lambdaGN1+1/3*lambdacN1*CN1, 'strain')) * detJ * gaussWeight(k);
        Rc = Rc + M(k, :)' * (DW_c - lambdacN1) * detJ * gaussWeight(k);
        RLambdaCv = RLambdaCv + kron(M(k, :)', matrixToVoigt(CxN1-CN1, 'strain')) * detJ * gaussWeight(k);
        RLambdaGv = RLambdaGv + kron(M(k, :)', matrixToVoigt(0.5*wedge(CN1, CN1)-GN1, 'strain')) * detJ * gaussWeight(k);
        RLambdac = RLambdac + M(k, :)' * (1 / 3 * innerProduct(GN1, CN1) - cN1) * detJ * gaussWeight(k);

        %% Tangent
        % KXX
        geometricalTangentPart = 2 * dNx' * lambdaCN1 * dNx;
        geometricalTangent = zeros(numberOfXDOFs);
        for g = 1:dimension
            geometricalTangent(g:dimension:numberOfXDOFs, g:dimension:numberOfXDOFs) = geometricalTangentPart;
        end
        KXX = KXX + geometricalTangent * detJ * gaussWeight(k);

        % KXLambdaC
        KXLambdaC = KXLambdaC + kron(M(k, :), 2*BN1') * detJ * gaussWeight(k);

        % =====

        % KTT
        dEtaN1_thetaN1 = kappa / thetaN1e;
        kTT1 = N(k, :)' * (etaN1 + thetaN1e * dEtaN1_thetaN1) / DT * N(k, :);
        kTT2 = -N(k, :)' * (etaN1 + (thetaN1e-thetaNe) * dEtaN1_thetaN1) / DT * N(k, :);
        kTT3 = dNx' * k0 * cN1^(-1) * GN1 * dNx;
        KTT = KTT + (kTT1 + kTT2 + kTT3) * detJ * gaussWeight(k);

        % KTG
        KTG = KTG + kron(M(k, :), dNx'*k0*1/cN1*BErs(dNx*thetaN1')) * detJ * gaussWeight(k);

        % KTc
        dEtaN1_cN1 = dimension * beta * c2;
        kTc1 = kron(N(k, :)'*M(k, :), (thetaN1e * dEtaN1_cN1 - (thetaN1e-thetaNe) * dEtaN1_cN1)/DT);
        kTc2 = -kron(M(k, :), dNx'*k0*1/(cN1^2)*GN1*(dNx * thetaN1'));
        KTc = KTc + (kTc1 + kTc2) * detJ * gaussWeight(k);

        % =====

        % KPD
        KPD = KPD + kron(M(k, :), dNx') * detJ * gaussWeight(k);

        % =====

        % KDP
        KDP = KDP + kron(M(k, :)', dNx) * detJ * gaussWeight(k);

        % KDD
        KDD = KDD + kron(M(k, :)'*M(k, :), 1/(er * e0 * sqrt(cN1))*CN1) * detJ * gaussWeight(k);

        % KDC
        KDC = KDC + kron(M(k, :)'*M(k, :), 1/(er * e0 * sqrt(cN1))*BErs(DN1)) * detJ * gaussWeight(k);

        % KDc
        D2W_D_cN1 = -1 / (2 * er * e0 * (cN1)^(3 / 2)) * CN1 * DN1;
        KDc = KDc + kron(M(k, :)'*M(k, :), D2W_D_cN1) * detJ * gaussWeight(k);

        % =====

        % KCD
        KCD = KCD + kron(M(k, :)'*M(k, :), 1/(2 * er * e0 * sqrt(cN1))*AErs(DN1)) * detJ * gaussWeight(k);

        % KCC
        KCC = KCC + kron(M(k, :)'*M(k, :), secDiffOperator(lambdaGN1)) * detJ * gaussWeight(k);

        % KCG
        KCG = KCG + 1 / 3 * kron(M(k, :)'*M(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

        % KCc
        D2W_C_cN1 = -1 / (4 * er * e0 * cN1^(3 / 2)) * matrixToVoigt(DN1*DN1', 'strain');
        KCc = KCc + kron(M(k, :)'*M(k, :), D2W_C_cN1) * detJ * gaussWeight(k);

        % KCLambdaC
        KCLambdaC = KCLambdaC - kron(M(k, :)'*M(k, :), map.Isym) * detJ * gaussWeight(k);

        % KCLambdaG
        KCLambdaG = KCLambdaG + kron(M(k, :)'*M(k, :), secDiffOperator(CN1)) * detJ * gaussWeight(k);

        % KCLambdac
        KCLambdac = KCLambdac + 1 / 3 * kron(M(k, :)'*M(k, :), matrixToVoigt(GN1, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KGC
        KGC = KGC + 1 / 3 * kron(M(k, :)'*M(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

        % KGLambdaG
        KGLambdaG = KGLambdaG - kron(M(k, :)'*M(k, :), map.Isym) * detJ * gaussWeight(k);

        % KGLambdac
        KGLambdac = KGLambdac + 1 / 3 * kron(M(k, :)'*M(k, :), matrixToVoigt(CN1, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KcD
        KcD = KcD - kron(M(k, :)'*M(k, :), 1/(2 * er * e0 * cN1^(3 / 2))*(CN1 * DN1)') * detJ * gaussWeight(k);

        % KcT
        D2W_c_thetaN1 = -dimension * beta * c2;
        KcT = KcT + kron(M(k, :)'*N(k, :), D2W_c_thetaN1) * detJ * gaussWeight(k);

        % KcC
        KcC = KcC - kron(M(k, :)'*M(k, :), 1/(4 * er * e0 * cN1^(3 / 2))*DN1'*BErs(DN1)*0.5) * detJ * gaussWeight(k);

        % Kcc
        D2W_c_cN1 = d1 / (2 * cN1^2) + c1 / (4 * cN1^(3 / 2)) + 3 / (8 * er * e0 * cN1^(5 / 2)) * DN1' * (CN1 * DN1);
        Kcc = Kcc + kron(M(k, :)'*M(k, :), D2W_c_cN1) * detJ * gaussWeight(k);

        % KcLambdac
        KcLambdac = KcLambdac - M(k, :)' * M(k, :) * detJ * gaussWeight(k);

        % =====

        % KLambdaCX
        KLambdaCX = KLambdaCX + 2 * kron(M(k, :)', BN1) * detJ * gaussWeight(k);

        % KLambdaCC
        KLambdaCC = KLambdaCC - kron(M(k, :)'*M(k, :), map.Isym) * detJ * gaussWeight(k);

        % =====

        % KLambdaGC
        KLambdaGC = KLambdaGC + kron(M(k, :)'*M(k, :), secDiffOperator(CN1)) * detJ * gaussWeight(k);

        % KLambdaGG
        KLambdaGG = KLambdaGG - kron(M(k, :)'*M(k, :), map.Isym) * detJ * gaussWeight(k);

        % =====

        % KLambdacC
        KLambdacC = KLambdacC + 1 / 3 * kron(M(k, :)'*M(k, :), matrixToVoigt(GN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacG
        KLambdacG = KLambdacG + 1 / 3 * kron(M(k, :)'*M(k, :), matrixToVoigt(CN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacc
        KLambdacc = KLambdacc - M(k, :)' * M(k, :) * detJ * gaussWeight(k);
    else

        %% Stress computation
        SN1 = 2 * lambdaCN1;
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        stressTensor.D = DN1;
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
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

function D = secDiffOperator(A)
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