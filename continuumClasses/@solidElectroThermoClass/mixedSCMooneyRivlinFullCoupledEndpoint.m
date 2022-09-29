function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinFullCoupledEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% MIXEDSCMOONEYRIVLINFULLCOUPLEDENDPOINT Element routine of class solidElectroThermoClass.
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
% In contrast to routine mixedSCMooneyRivlinEndpoint, this routine
% also considers thermoelectric coupling.
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
%        every field, here: (X, T, D, C, G, c, LambdaC, LambdaG, Lambdac)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, T, D, C, G, c, LambdaC,
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
% In preparation
%
% SEE ALSO
% mixedSCMooneyRivlinMidpoint
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
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
if strcmpi(mixedFEObject.typeShapeFunction,'sameOrder')
    M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
elseif strcmpi(mixedFEObject.typeShapeFunction,'detailedOrder')
    M_k_I = mixedFEObject.shapeFunctionObject.N_k_I{1};
end
% material data and voigt notation
numberOfInternalNodes = size(M_k_I, 2);
dimension = obj.dimension;
DT = setupObject.timeStepSize;
time = setupObject.time;

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
edN = obj.qN(edof, 1:dimension).';
phiN = obj.qN(edof, dimension+1).';
edAlphaN = obj.mixedFEObject.qN(e, :);
thetaN = obj.qN(edof, dimension+2).';

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
rho0 = materialObject.rhoSource;
R = materialObject.RSource;

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
    [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % assembly internal dofs
    DN = reshape(extractedDN, 3, []) * M_k_I(k, :)';
    CNv = reshape(extractedCNv, 6, []) * M_k_I(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    GNv = reshape(extractedGNv, 6, []) * M_k_I(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    cN = extractedcN.' * M_k_I(k, :)';
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
    % Temperature
    thetaN1e = N_k_I(k, :) * thetaN1.';
    thetaNe = N_k_I(k, :) * thetaN.';
    % Deformation Gradient
    FxN1 = edN1 * dN_X_I';
    CxN1 = FxN1.' * FxN1;
    % Electrical quantities
    EN1 = -dN_X_I * phiN1.';
    EN = -dN_X_I * phiN.';
    % B-matrix (endpoint configuration)
    BN1 = BMatrix(dN_X_I, FxN1);

    if ~computePostData

        %% energy
        fThetaN = thetaNe / thetaR;
        fThetaN1 = thetaN1e / thetaR;
        DfThetaN_theta = 1 / thetaR;
        DfThetaN1_theta = 1 / thetaR;

        WTN = kappa * (thetaNe - thetaR - thetaNe * log(thetaNe/thetaR));
        WTMN = -dimension * beta * c2 * (thetaNe - thetaR) * (cN - 1);
        WEMN = 1 / (2 * er * e0 * cN^(1 / 2)) * DN' * CN * DN;
        WMIsoN = a * (trace(CN) - 3) + b * (trace(GN) - 3);
        WMVolN = c1 / 2 * (sqrt(cN) - 1)^2 - d1 * log(sqrt(cN));
        WN = WTN + WTMN + fThetaN * (WEMN + WMIsoN + WMVolN);
        DWTN_thetaN = -kappa * log(thetaNe/thetaR);
        DWTMN_thetaN = -dimension * beta * c2 * (cN - 1);
        DW_thetaN = DWTN_thetaN + DWTMN_thetaN + DfThetaN_theta * (WEMN + WMIsoN + WMVolN);
        etaN = -DW_thetaN;
        u0N = WN - DN' * EN + thetaNe * etaN;

        WTN1 = kappa * (thetaN1e-thetaR - thetaN1e * log(thetaN1e/thetaR));
        WTMN1 = -dimension * beta * c2 * (thetaN1e-thetaR) * (cN1 - 1);
        WEMN1 = 1 / (2 * er * e0 * cN1^(1 / 2)) * DN1.' * CN1 * DN1;
        WMIsoN1 = a * (trace(CN1) - 3) + b * (trace(GN1) - 3);
        WMVolN1 = c1 / 2 * (sqrt(cN1) - 1)^2 - d1 * log(sqrt(cN1));
        WN1 = WTN1 + WTMN1 + fThetaN1 * (WEMN1 + WMIsoN1 + WMVolN1);
        DWTN1_theta = -kappa * log(thetaN1e/thetaR);
        DWTMN1_theta = -dimension * beta * c2 * (cN1 - 1);
        DW_theta = DWTN1_theta + DWTMN1_theta + DfThetaN1_theta * (WEMN1 + WMIsoN1 + WMVolN1);
        etaN1 = -DW_theta;
        u0N1 = WN1 - DN1.' * EN1 + thetaN1e * etaN1;

        elementEnergy.internalEnergy = elementEnergy.internalEnergy + u0N1 * detJ * gaussWeight(k);
        elementEnergy.helmholtz = elementEnergy.helmholtz + WN1 * detJ * gaussWeight(k);
        elementEnergy.DE = elementEnergy.DE + DN1.' * EN1 * detJ * gaussWeight(k);
        elementEnergy.TS = elementEnergy.TS + thetaN1e * etaN1 * detJ * gaussWeight(k);
        elementEnergy.S = elementEnergy.S + etaN1 * detJ * gaussWeight(k);
        elementEnergy.deltaS = elementEnergy.deltaS + (etaN1 - etaN) * detJ * gaussWeight(k);
        elementEnergy.deltaU = elementEnergy.deltaU + (u0N1 - u0N) * detJ * gaussWeight(k);

        %% Derivatives
        DW_theta = DfThetaN1_theta * (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) + c1 / 2 * (sqrt(cN1) - 1)^2 - d1 * log(sqrt(cN1)) + 1 / (2 * er * e0 * sqrt(cN1)) * DN1.' * (CN1 * DN1)) - 3 * beta * c2 * (cN1 - 1) - kappa * log(thetaN1e/thetaR);
        DW_D = fThetaN1 * 1 / (er * e0 * sqrt(cN1)) * CN1 * DN1;
        DW_C = fThetaN1 * (a * eye(3) + 1 / (2 * er * e0 * sqrt(cN1)) * (DN1 * DN1.'));
        DW_G = fThetaN1 * b * eye(3);
        DW_c = fThetaN1 * (-d1 / (2 * cN1) + c1 / 2 * (1 - 1 / sqrt(cN1)) - 1 / (4 * er * e0 * cN1^(3 / 2)) * DN1.' * (CN1 * DN1)) - dimension * beta * c2 * (thetaN1e-thetaR);

        %% Constitutive equations
        etaN1 = -DW_theta;
        QN1 = -k0 * cN1^(-1) * GN1 * (dN_X_I * thetaN1.');

        %% Sources
        if isa(rho0,'function_handle')
            XGP = kron(N_k_I(k,:),eye(dimension))*edR(:);
            rho0GP = rho0(XGP(1),XGP(2),XGP(3));
        else
            rho0GP = rho0;
        end    
        if isa(R,'function_handle')
            XGP = kron(N_k_I(k,:),eye(dimension))*edR(:);
            RGP = R(XGP(1),XGP(2),XGP(3));
            flagTransient = false;
        else
            flagTransient = true;
            RGP = R;
        end

        %% Residual
        RX = RX + BN1.' * 2 * matrixToVoigt(lambdaCN1, 'stress') * detJ * gaussWeight(k);
        if flagTransient
            RT = RT + (N_k_I(k, :)' * (thetaN1e * etaN1 - thetaNe * etaN) / DT + N_k_I(k, :)' * (thetaN1e-thetaNe) / DT * DW_theta - dN_X_I' * QN1 - N_k_I(k,:)'*RGP) * detJ * gaussWeight(k);
        else
            RT = RT + (-dN_X_I' * QN1 - N_k_I(k,:)'*materialObject.timeFunctionRSource(time)*RGP) * detJ * gaussWeight(k);
        end
        RP = RP + (dN_X_I' * DN1 + materialObject.timeFunctionRhoSource(time)*N_k_I(k,:)'*rho0GP) * detJ * gaussWeight(k);
        RD = RD + kron(M_k_I(k, :)', (DW_D + dN_X_I * phiN1.')) * detJ * gaussWeight(k);
        RCv = RCv + kron(M_k_I(k, :)', matrixToVoigt(DW_C-lambdaCN1+wedge(lambdaGN1, CN1)+1/3*lambdacN1*GN1, 'strain')) * detJ * gaussWeight(k);
        RGv = RGv + kron(M_k_I(k, :)', matrixToVoigt(DW_G-lambdaGN1+1/3*lambdacN1*CN1, 'strain')) * detJ * gaussWeight(k);
        Rc = Rc + M_k_I(k, :)' * (DW_c - lambdacN1) * detJ * gaussWeight(k);
        RLambdaCv = RLambdaCv + kron(M_k_I(k, :)', matrixToVoigt(CxN1-CN1, 'strain')) * detJ * gaussWeight(k);
        RLambdaGv = RLambdaGv + kron(M_k_I(k, :)', matrixToVoigt(0.5*wedge(CN1, CN1)-GN1, 'strain')) * detJ * gaussWeight(k);
        RLambdac = RLambdac + M_k_I(k, :)' * (1 / 3 * innerProduct(GN1, CN1) - cN1) * detJ * gaussWeight(k);

        %% Tangent
        % KXX
        geometricalTangentPart = 2 * dN_X_I' * lambdaCN1 * dN_X_I;
        geometricalTangent = zeros(numberOfXDOFs);
        for g = 1:dimension
            geometricalTangent(g:dimension:numberOfXDOFs, g:dimension:numberOfXDOFs) = geometricalTangentPart;
        end
        KXX = KXX + geometricalTangent * detJ * gaussWeight(k);

        % KXLambdaC
        KXLambdaC = KXLambdaC + kron(M_k_I(k, :), 2*BN1') * detJ * gaussWeight(k);

        % =====

        % KTT
        D2W_theta_thetaN1 = -kappa / thetaN1e;
        dEtaN1_thetaN1 = kappa / thetaN1e;
        if flagTransient
            kTT1 = N_k_I(k, :)' * (etaN1 + thetaN1e * dEtaN1_thetaN1) / DT * N_k_I(k, :);
            kTT2 = N_k_I(k, :)' * (DW_theta + (thetaN1e-thetaNe) * D2W_theta_thetaN1) / DT * N_k_I(k, :);
        else 
            kTT1 = 0;
            kTT2 = 0;
        end
        kTT3 = dN_X_I' * k0 * cN1^(-1) * GN1 * dN_X_I;
        KTT = KTT + (kTT1 + kTT2 + kTT3) * detJ * gaussWeight(k);

        % KTD
        if flagTransient
            kTD1 = -kron(N_k_I(k, :)'*M_k_I(k, :), 1/DT*thetaN1e*DfThetaN1_theta*1/(er * e0 * sqrt(cN1))*(CN1 * DN1)');
            kTD2 = kron(N_k_I(k, :)'*M_k_I(k, :), DfThetaN1_theta*1/DT*(thetaN1e-thetaNe)*1/(er * e0 * sqrt(cN1))*(CN1 * DN1)');
        else
            kTD1 = 0;
            kTD2 = 0;
        end
        KTD = KTD + (kTD1 + kTD2) * detJ * gaussWeight(k);

        % KTC
        D2W_theta_CN1 = DfThetaN1_theta * (matrixToVoigt(a*eye(3), 'strain')' + 1 / (2 * er * e0 * sqrt(cN1)) * DN1' * BErs(DN1));
        if flagTransient
            kTC1 = kron(N_k_I(k, :)'*M_k_I(k, :), (-thetaN1e * D2W_theta_CN1 + (thetaN1e-thetaNe) * D2W_theta_CN1)/DT);
        else
            kTC1 = 0;
        end
        KTC = KTC + (kTC1) * detJ * gaussWeight(k);

        % KTG
        D2W_theta_GN1 = matrixToVoigt(DfThetaN1_theta*b*eye(3), 'strain')';
        if flagTransient
            kTG1 = kron(N_k_I(k, :)'*M_k_I(k, :), (-thetaN1e * D2W_theta_GN1 + (thetaN1e-thetaNe) * D2W_theta_GN1)/DT);
        else
            kTG1 = 0;
        end
        kTG2 = kron(M_k_I(k, :), dN_X_I'*k0*1/cN1*BErs(dN_X_I*thetaN1'));
        KTG = KTG + (kTG1 + kTG2) * detJ * gaussWeight(k);

        % KTc
        D2W_theta_cN1 = DfThetaN1_theta * (c1 / 2 * (1 - 1 / sqrt(cN1)) - d1 / (2 * cN1) - 1 / (4 * er * e0 * (cN1)^(3 / 2)) * DN1' * (CN1 * DN1)) - dimension * beta * c2;
        if flagTransient
            kTc1 = kron(N_k_I(k, :)'*M_k_I(k, :), (-thetaN1e * D2W_theta_cN1 + (thetaN1e-thetaNe) * D2W_theta_cN1)/DT);
        else
            kTc1 = 0;
        end
        kTc2 = -kron(M_k_I(k, :), dN_X_I'*k0*1/(cN1^2)*GN1*(dN_X_I * thetaN1'));
        KTc = KTc + (kTc1 + kTc2) * detJ * gaussWeight(k);

        % =====

        % KPD
        KPD = KPD + kron(M_k_I(k, :), dN_X_I') * detJ * gaussWeight(k);

        % =====

        % KDT
        KDT = KDT + kron(M_k_I(k, :)'*N_k_I(k, :), DfThetaN1_theta*1/(er * e0 * sqrt(cN1))*CN1*DN1) * detJ * gaussWeight(k);

        % KDP
        KDP = KDP + kron(M_k_I(k, :)', dN_X_I) * detJ * gaussWeight(k);

        % KDD
        KDD = KDD + kron(M_k_I(k, :)'*M_k_I(k, :), fThetaN1*1/(er * e0 * sqrt(cN1))*CN1) * detJ * gaussWeight(k);

        % KDC
        KDC = KDC + kron(M_k_I(k, :)'*M_k_I(k, :), fThetaN1*1/(er * e0 * sqrt(cN1))*BErs(DN1)) * detJ * gaussWeight(k);

        % KDc
        D2W_D_cN1 = -fThetaN1 * 1 / (2 * er * e0 * (cN1)^(3 / 2)) * CN1 * DN1;
        KDc = KDc + kron(M_k_I(k, :)'*M_k_I(k, :), D2W_D_cN1) * detJ * gaussWeight(k);

        % =====

        % KCT
        KCT = KCT + kron(M_k_I(k, :)'*N_k_I(k, :), DfThetaN1_theta*matrixToVoigt(a*eye(3)+1/(2 * er * e0 * sqrt(cN1))*(DN1 * DN1'), 'strain')) * detJ * gaussWeight(k);

        % KCD
        KCD = KCD + kron(M_k_I(k, :)'*M_k_I(k, :), fThetaN1*1/(2 * er * e0 * sqrt(cN1))*AErs(DN1)) * detJ * gaussWeight(k);

        % KCC
        KCC = KCC + kron(M_k_I(k, :)'*M_k_I(k, :), secDiffOperator(lambdaGN1)) * detJ * gaussWeight(k);

        % KCG
        KCG = KCG + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

        % KCc
        D2W_C_cN1 = - fThetaN1 / (4 * er * e0 * cN1^(3 / 2)) * matrixToVoigt(DN1*DN1', 'strain');
        KCc = KCc + kron(M_k_I(k, :)'*M_k_I(k, :), D2W_C_cN1) * detJ * gaussWeight(k);

        % KCLambdaC
        KCLambdaC = KCLambdaC - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

        % KCLambdaG
        KCLambdaG = KCLambdaG + kron(M_k_I(k, :)'*M_k_I(k, :), secDiffOperator(CN1)) * detJ * gaussWeight(k);

        % KCLambdac
        KCLambdac = KCLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN1, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KGT
        KGT = KGT + kron(M_k_I(k, :)'*N_k_I(k, :), DfThetaN1_theta*b*matrixToVoigt(eye(3), 'strain')) * detJ * gaussWeight(k);

        % KGC
        KGC = KGC + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), lambdacN1*map.Isym) * detJ * gaussWeight(k);

        % KGLambdaG
        KGLambdaG = KGLambdaG - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

        % KGLambdac
        KGLambdac = KGLambdac + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN1, 'strain')) * detJ * gaussWeight(k);

        % =====

        % KcT
        D2W_c_thetaN1 = DfThetaN1_theta * (c1 / 2 * (1 - 1 / sqrt(cN1)) - d1 / (2 * cN1) - 1 / (4 * er * e0 * (cN1)^(3 / 2)) * DN1' * (CN1 * DN1)) - dimension * beta * c2;
        KcT = KcT + kron(M_k_I(k, :)'*N_k_I(k, :), D2W_c_thetaN1) * detJ * gaussWeight(k);

        % KcD
        KcD = KcD - kron(M_k_I(k, :)'*M_k_I(k, :), fThetaN1*1/(2 * er * e0 * cN1^(3 / 2))*(CN1 * DN1)') * detJ * gaussWeight(k);

        % KcC
        KcC = KcC - kron(M_k_I(k, :)'*M_k_I(k, :), fThetaN1*1/(4 * er * e0 * cN1^(3 / 2))*DN1'*BErs(DN1)) * detJ * gaussWeight(k);

        % Kcc
        D2W_c_cN1 = fThetaN1 * (d1 / (2 * cN1^2) + c1 / (4 * cN1^(3 / 2)) + 3 / (8 * er * e0 * cN1^(5 / 2)) * DN1' * (CN1 * DN1));
        Kcc = Kcc + kron(M_k_I(k, :)'*M_k_I(k, :), D2W_c_cN1) * detJ * gaussWeight(k);

        % KcLambdac
        KcLambdac = KcLambdac - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);

        % =====

        % KLambdaCX
        KLambdaCX = KLambdaCX + 2 * kron(M_k_I(k, :)', BN1) * detJ * gaussWeight(k);

        % KLambdaCC
        KLambdaCC = KLambdaCC - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

        % =====

        % KLambdaGC
        KLambdaGC = KLambdaGC + kron(M_k_I(k, :)'*M_k_I(k, :), secDiffOperator(CN1)) * detJ * gaussWeight(k);

        % KLambdaGG
        KLambdaGG = KLambdaGG - kron(M_k_I(k, :)'*M_k_I(k, :), map.Isym) * detJ * gaussWeight(k);

        % =====

        % KLambdacC
        KLambdacC = KLambdacC + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(GN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacG
        KLambdacG = KLambdacG + 1 / 3 * kron(M_k_I(k, :)'*M_k_I(k, :), matrixToVoigt(CN1, 'strain')') * detJ * gaussWeight(k);

        % KLambdacc
        KLambdacc = KLambdacc - M_k_I(k, :)' * M_k_I(k, :) * detJ * gaussWeight(k);
    else

        %% Stress computation
        SN1 = 2 * lambdaCN1;
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        stressTensor.D = DN1;
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
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