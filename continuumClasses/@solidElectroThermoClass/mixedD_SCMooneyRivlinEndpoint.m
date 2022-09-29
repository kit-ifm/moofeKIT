function [rData, kData, elementEnergy, array] = mixedD_SCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
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
edof = meshObject.edof(e,:);
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
if strcmpi(mixedFEObject.typeShapeFunction,'sameOrder')
    N_D = mixedFEObject.shapeFunctionObject.N;
elseif strcmpi(mixedFEObject.typeShapeFunction,'detailedOrder')
    N_D = mixedFEObject.shapeFunctionObject.N{1};
end
numberOfInternalNodes = size(N_D, 2);
% material data and voigt notation
dimension = obj.dimension;
% 
edN1 = dofs.edN1;
phiN1 = dofs.phiN1;
edAlphaN1 = dofs.edAlphaN1;
thetaN1 = dofs.thetaN1;
% 
edN  = obj.qN(edof,1:dimension)';
phiN = obj.qN(edof,dimension+1)';
edAlphaN = obj.mixedFEObject.qN(e, :);
thetaN = obj.qN(edof,dimension+2)';
% 
extractedDN = edAlphaN(1:3*numberOfInternalNodes).';
extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';

edRef = obj.qR(edof,1:dimension)';
%% Jacobi-Matrix
J = edRef*dNr';
JN1 = edN1*dNr';
%% Energy
elementEnergy.internalEnergy = 0;
elementEnergy.helmholtz = 0;
elementEnergy.DE = 0;
elementEnergy.TS = 0;
elementEnergy.S = 0;
elementEnergy.deltaS = 0;
elementEnergy.deltaU = 0;
%% Material parameters
a = materialObject.a;
b = materialObject.b;
c1 = materialObject.c1;
d1 = materialObject.d1;
c2 = materialObject.c2;
d2 = materialObject.d2;
e0 = materialObject.e0;
er = materialObject.e1;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;
rho0 = materialObject.rhoSource;
R = materialObject.RSource;
DT = setupObject.timeStepSize;
time = setupObject.time;
% map
map.Voigt = [1, 5, 9, 4, 8, 7]';
map.VoigtInv = [1, 4, 6; 4, 2, 5; 6, 5, 3];
map.VoigtFull = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
map.VoigtFullInv = [1, 4, 6; 7, 2, 5; 9, 8, 3];
map.Isym = [eye(3), zeros(3); zeros(3), 2 * eye(3)];
map.Isyminv = [eye(3), zeros(3); zeros(3), 1 / 2 * eye(3)];
%% Initialization
RX = rData{1,1};
RP = rData{2,1};
RT = rData{3,1};
RD = rData{4,1};
KXX = kData{1,1}; KXP = kData{1,2}; KXT = kData{1,3}; KXD = kData{1,4};
KPX = kData{2,1}; KPP = kData{2,2}; KPT = kData{2,3}; KPD = kData{2,4};
KTX = kData{3,1}; KTP = kData{3,2}; KTT = kData{3,3}; KTD = kData{3,4};
KDX = kData{4,1}; KDP = kData{4,2}; KDT = kData{4,3}; KDD = kData{4,4};

%% DOFs
numberOfMechanicalDOFs = numel(RX);
numberOfElectricalDOFs = numel(RP);
numberOfInternalDOFs = numel(RD);
numberOfTemperatureDOFs = numel(RT);
numberOfInternalNodes = size(obj.mixedFEObject.shapeFunctionObject.N,2);
%% Gauss loop
for k = 1:numberOfGausspoints
    index = dimension*k-(dimension-1):dimension*k;
    detJ = det(J(:,index)');
    detJN1 = det(JN1(:,index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,index)')\dNr(index,:);
    % Temperature
    thetaN1e = N(k,:)*thetaN1.';
    thetaNe = N(k,:)*thetaN.';
    DotTheta = (thetaN1e-thetaNe)/DT;  
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
    CxN1 = FxN1.'*FxN1;
    CxN = FxN.'*FxN;
    DotCx = (CxN1-CxN)/DT;
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN = 0.5*wedge(CxN,CxN);
    cxN1 = det(CxN1);
    cxN = det(CxN);
    % Electrical quantities
    EN1 = -dNx*phiN1.';
    EN = -dNx*phiN.';
%     DN1 = zeros(3,1);
%     for m = 1:numberOfInternalNodes
%         DN1 = DN1 + N_D(k,m)*DN1e(e,ones(1,3)*(m - 1)*3 + (1:3))';
%     end
    DN1 = reshape(extractedDN1, 3, []) * N_D(k, :)';
%     DN = zeros(3,1);
%     for m = 1:numberOfInternalNodes
%         DN = DN + N_D(k,m)*DNe(e,ones(1,3)*(m - 1)*3 + (1:3))';
%     end
    DN = reshape(extractedDN, 3, []) * N_D(k, :)';
    % B-matrix (current configuration)
    BN1 = BMatrix(dNx,FxN1);
    BN1 = 2*BN1;
    %% energy
    PsiTN = kappa*(thetaNe - thetaR - thetaNe*log(thetaNe/thetaR));
    PsiTMN = -dimension*beta*c2*(thetaNe - thetaR)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiMIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiMVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
    s0N = mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaNe,3);
    u0N = PsiN - DN'*EN + thetaNe*s0N;
    PsiTN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
%     PsiTMN1 = DIM*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiTMN1 = -dimension*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiMIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiMVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
    s0N1 = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,3);
    u0N1 = PsiN1 - DN1'*EN1 + thetaN1e*s0N1;
    % MR-Mike
%     Psi = a*(trace(CN1)-DIM) + b*(1/2*((trace(CN1))^2 - trace((CN1)^2))-DIM) - d1*log(sqrt(cxN1)) + kappa*(thetaN1 - thetaR -thetaN1*log(thetaN1/thetaR)) + DIM*beta*(thetaN1 - thetaR)*d1*I3N1^(-1);
    % MR-ThermoMech
%     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR)) - DIM*beta*(thetaN1e - thetaR)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    % MR-ThermoElectroMech
%     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR)) - DIM*beta*(thetaN1e - thetaR)*c2*(cxN1 - 1) + 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
%     
    elementEnergy.internalEnergy = elementEnergy.internalEnergy + u0N1*detJ*gaussWeight(k);
    elementEnergy.helmholtz = elementEnergy.helmholtz + PsiN1*detJ*gaussWeight(k);
    elementEnergy.DE = elementEnergy.DE + DN1'*EN1*detJ*gaussWeight(k);
    elementEnergy.TS = elementEnergy.TS + thetaN1e*s0N1*detJ*gaussWeight(k);
    elementEnergy.S = elementEnergy.S + s0N1*detJ*gaussWeight(k);
    elementEnergy.deltaS = elementEnergy.deltaS + (s0N1-s0N)*detJ*gaussWeight(k);
    elementEnergy.deltaU = elementEnergy.deltaU + (u0N1-u0N)*detJ*gaussWeight(k);
    %% residual
    DPsi_CN1 = 1/(2*er*e0*cxN1^(1/2))*DN1*DN1.' + a*eye(3);
    DPsi_GN1 = b*eye(3);
    DPsi_cN1 = -dimension*beta*c2*(thetaN1e - thetaR) - 1/(4*er*e0*cxN1^(3/2))*DN1.'*CxN1*DN1 + c1/2*(1-cxN1^(-1/2)) - d1/(2*cxN1);
    SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
    if ~computePostData
        QxN1 = -k0*1/cxN1*GxN1*(dNx*thetaN1.');
        DPsi_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;
        %% Sources
        if isa(rho0,'function_handle')
            XGP = kron(N(k,:),eye(dimension))*edRef(:);
            rho0GP = rho0(XGP(1),XGP(2),XGP(3));
        else
            rho0GP = rho0;
        end    
        if isa(R,'function_handle')
            XGP = kron(N(k,:),eye(dimension))*edRef(:);
            RGP = R(XGP(1),XGP(2),XGP(3));
            flagTransient = false;
        else
            flagTransient = true;
            RGP = R;
        end

        %% Residual
        RX = RX + BN1.'*0.5*SN1(map.Voigt)*detJ*gaussWeight(k);
        RP = RP + (dNx'*DN1 + materialObject.timeFunctionRhoSource(time)*N(k,:)'*rho0GP)*detJ*gaussWeight(k);
        if flagTransient
            RT = RT + (N(k,:)'*(thetaN1e*s0N1-thetaNe*s0N)/DT - N(k,:)'*(thetaN1e - thetaNe)/DT*s0N1 - dNx'*QxN1 - N(k,:)'*RGP)*detJ*gaussWeight(k);
        else
            RT = RT + (-dNx'*QxN1 - materialObject.timeFunctionRSource(time)*N(k,:)'*RGP)*detJ*gaussWeight(k);
        end
        RD = RD + kron(N_D(k,:)',(DPsi_DN1 - EN1))*detJ*gaussWeight(k);        
        
        %% Tangent
        % ---------------------------- KXX ----------------------------
        D2Psi_CN1cN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        D2Psi_cN1cN1 = 3/(8*er*e0)*cxN1^(-5/2)*DN1'*(CxN1*DN1) + c1/(4*cxN1^(3/2)) + d1/(2*cxN1^2);
        D2Psi_cN1CN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        elasticityTensor = 4*D2Psi_CN1cN1(map.Voigt)*GxN1(map.Voigt)' + 2*secondDiffFunction(DPsi_GN1) + 2*DPsi_cN1*secondDiffFunction(CxN1) + 4*D2Psi_cN1cN1*(GxN1(map.Voigt)*GxN1(map.Voigt)') + 4*GxN1(map.Voigt)*D2Psi_cN1CN1(map.Voigt)';
        materialTangent = 0.25*BN1'*elasticityTensor*BN1;
        % Geometrical part of tangent operator
        geometricalTangentPart = dNx'*SN1*dNx;
        geometricalTangent = zeros(numberOfMechanicalDOFs);
        for g = 1:dimension
            geometricalTangent(g:dimension:numberOfMechanicalDOFs,g:dimension:numberOfMechanicalDOFs) = geometricalTangentPart;
        end
        KXX = KXX + (materialTangent + geometricalTangent)*detJ*gaussWeight(k);
        % ---------------------------- KXP ----------------------------
        KXP = zeros(numberOfMechanicalDOFs,numberOfElectricalDOFs);
        % ---------------------------- KXT ----------------------------
        D2Psi_cN1thetaN1 = -dimension*beta*c2;
        DSN1_thetaN1 = 2*D2Psi_cN1thetaN1*GxN1;
        KXT = KXT + (BN1'*0.5*DSN1_thetaN1(map.Voigt))*N(k,:)*detJ*gaussWeight(k);
        % ---------------------------- KXD ----------------------------
        BEN1 = zeros(6,numberOfInternalDOFs);
        BEN1(1,1:3:end) = 2*DN1(1)*N_D(k,:);
        BEN1(2,2:3:end) = 2*DN1(2)*N_D(k,:);
        BEN1(3,3:3:end) = 2*DN1(3)*N_D(k,:);
        BEN1(4,1:3:end) = DN1(2)*N_D(k,:);
        BEN1(4,2:3:end) = DN1(1)*N_D(k,:);
        BEN1(5,2:3:end) = DN1(3)*N_D(k,:);
        BEN1(5,3:3:end) = DN1(2)*N_D(k,:);
        BEN1(6,1:3:end) = DN1(3)*N_D(k,:);
        BEN1(6,3:3:end) = DN1(1)*N_D(k,:);
        D2Psi_cN1DN1 = -1/(2*er*e0)*cxN1^(-3/2)*(CxN1*DN1);
        temp1 = BN1'*GxN1(map.Voigt);
        temp2 = kron(N_D(k,:),D2Psi_cN1DN1);
        KXD = KXD + (BN1'*0.5*2/(2*er*e0*cxN1^(1/2))*BEN1 + temp1*temp2(:)')*detJ*gaussWeight(k);
        % ---------------------------- KPX ----------------------------
        %
        % ---------------------------- KPP ----------------------------
        %
        % ---------------------------- KPT ----------------------------
        %
        % ---------------------------- KPD ----------------------------
        KPD = KPD + kron(N_D(k,:),dNx')*detJ*gaussWeight(k);
        % ---------------------------- KTX ----------------------------
        Ds0N1_cxN1 = dimension*beta*c2;
        if flagTransient
            kTX1 = N(k,:)'*thetaN1e*Ds0N1_cxN1/DT*(BN1'*GxN1(map.Voigt))';
            kTX2 = -N(k,:)'*(thetaN1e - thetaNe)*Ds0N1_cxN1/DT*(BN1'*GxN1(map.Voigt))';
        else
            kTX1 = 0;
            kTX2 = 0;
        end
        kTX3 = -dNx'*k0*cxN1^(-2)*(GxN1*(dNx*thetaN1'))*(BN1'*GxN1(map.Voigt))';
        kTX4 = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for jj = 1:numberOfMechanicalDOFs
            BN1jj = map.Isyminv*BN1(:,jj);
            kTX4(:,jj) = dNx'*k0*cxN1^(-1)*wedge(CxN1,BN1jj(map.VoigtInv))*(dNx*thetaN1');
        end
        KTX = KTX + (kTX1 + kTX2 + kTX3 + kTX4)*detJ*gaussWeight(k);
        % ---------------------------- KTP ----------------------------
        %
        % ---------------------------- KTT ----------------------------
        %
        % ---------------------------- KTD ----------------------------
        Ds0N1_thetaN1 = kappa/thetaN1e;
        DQN1_thetaN1 = -k0*(cxN1^(-1)*GxN1*dNx*eye(size(dNx,2)));
        if flagTransient
            kTT1 = N(k,:)'*N(k,:)*s0N1/DT - N(k,:)'*N(k,:)*s0N1/DT;
            kTT2 = N(k,:)'*(thetaN1e-(thetaN1e - thetaNe))*Ds0N1_thetaN1*N(k,:)/DT;
        else
            kTT1 = 0;
            kTT2 = 0;
        end
        kTT3 = - dNx'*DQN1_thetaN1;
        KTT = KTT + (kTT1 + kTT2 + kTT3)*detJ*gaussWeight(k);
        % ---------------------------- KDX ----------------------------
        KDX = KDX + (BEN1'*0.5*2*1/(2*er*e0*cxN1^(1/2))*BN1 + temp2(:)*temp1')*detJ*gaussWeight(k);
        % ---------------------------- KDP ----------------------------
        KDP = KDP + kron(N_D(k,:),dNx')'*detJ*gaussWeight(k);
        % ---------------------------- KDT ----------------------------
        %
        % ---------------------------- KDD ----------------------------
        D2Psi_DN1DN1 = 1/(er*e0*cxN1^(1/2))*CxN1;
        KDD = KDD + kron(N_D(k,:)'*N_D(k,:),D2Psi_DN1DN1)*detJ*gaussWeight(k);
    else
        %% Stress computation
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
    kData{1,1} = KXX;   kData{1,2} = KXP;   kData{1,3} = KXT;   kData{1,4} = KXD;
    kData{2,1} = KPX;   kData{2,2} = KPP;   kData{2,3} = KPT;   kData{2,4} = KPD;
    kData{3,1} = KTX;   kData{3,2} = KTP;   kData{3,3} = KTT;   kData{3,4} = KTD;
    kData{4,1} = KDX;   kData{4,2} = KDP;   kData{4,3} = KDT;   kData{4,4} = KDD;
end
end

function out = mooneyRivlin(obj,C,G,c,D,theta,choice)
a = obj.materialObject.a;
b = obj.materialObject.b;
c1 = obj.materialObject.c1;
d1 = obj.materialObject.d1;
c2 = obj.materialObject.c2;
d2 = obj.materialObject.d2;
e0 = obj.materialObject.e0;
er = obj.materialObject.e1;
kappa = obj.materialObject.kappa;
beta = obj.materialObject.beta;
thetaR = obj.materialObject.thetaR;
dimension = obj.dimension;
% 
PsiT = kappa*(theta - thetaR - theta*log(theta/thetaR));
PsiTM = -dimension*beta*c2*(theta - thetaR)*(c - 1);
PsiEM = 1/(2*er*e0*c^(1/2))*D'*C*D;
PsiMIso = a*(trace(C)-3) + b*(trace(G)-3);
PsiMVol = c1/2*(sqrt(c)-1)^2 - d1*log(sqrt(c));
Psi = PsiT + PsiTM + PsiEM + PsiMIso + PsiMVol;
% 
DPsi_C = 1/(2*er*e0)*c^(-1/2)*D*D' + a*eye(3);
DPsi_c = -dimension*beta*c2*(theta - thetaR) - 1/(4*er*e0*c^(3/2))*D'*C*D + c1/2*(1-c^(-1/2)) - d1/(2*c);
DPsi_D = 1/(er*e0*c^(1/2))*C*D;
DPsi_t = -kappa*log(theta/thetaR) - dimension*beta*c2*(c - 1);
% 
D2Psi_Cc = -1/(4*er*e0)*c^(-3/2)*D*D';
% 
s0 = -DPsi_t;
Ds0_c = dimension*beta*c2;
% 
% u0 = Psi + theta*s0;
% u0 = Psi + theta*s0 - D'*E;
u0 = kappa*(theta-thetaR) + dimension*beta*c2*thetaR*(c-1) + 1/(2*er*e0*c^(1/2))*D'*C*D + a*(trace(C)-3) + b*trace(G) - d1*log(c^(1/2)) + c1/2*(c^(1/2) - 1)^2;
Du0_C = 1/(2*er*e0*c^(1/2))*D*D' + a*eye(3);
Du0_c = dimension*beta*c2*thetaR - 1/(4*er*e0*c^(3/2))*D'*C*D + c1/2*(1 - c^(-1/2)) - d1/(2*c);
Du0_D = 1/(er*e0*c^(1/2))*C*D;
% 
if choice == 1
    out = Psi;
elseif choice == 2
    out = u0;
elseif choice == 3
    out = s0;
elseif choice == 4
    out = Du0_C;
elseif choice == 5
    out = Du0_c;
elseif choice == 6
    out = Du0_D;
elseif choice == 14
    out = DPsi_C;
elseif choice == 15
    out = DPsi_c;
elseif choice == 16
    out = DPsi_D;
elseif choice == 17
    out = DPsi_t;
elseif choice == 21
    out = D2Psi_Cc;
end
end

function [rData,kData,GlobalEnergy] = Residuum2(j,obj,map,time,DT,stress,indexD,rData,kData,dofsNT,flagNumericalTangent)
edofLocal = obj.EDOF{j};
%% Shape functions
N = obj.SHAPEF.N;
N_D = obj.SHAPEF_D.N;
wp = obj.SHAPEF.wp;
NGP = obj.NGP;
dNr = obj.SHAPEF.dNr;
DIM = obj.DIM;
% 
if flagNumericalTangent
    edN1 = dofsNT.edN1; 
    phiN1 = dofsNT.phiN1;
    DN1e = dofsNT.DN1;
    thetaN1 = dofsNT.thetaN1;
else
    edN1  = obj.QN1(edofLocal,1:DIM)';
    phiN1 = obj.QN1(edofLocal,DIM+1)';
    DN1e  = obj.DN1;
    thetaN1 = obj.QN1(edofLocal,DIM+2)';
end    
% 
edN  = obj.QN(edofLocal,1:DIM)';
phiN = obj.QN(edofLocal,DIM+1)';
DNe  = obj.DN;
thetaN = obj.QN(edofLocal,DIM+2)';
% 
edRef = obj.QREF(edofLocal,1:DIM)';
%% Jacobi-Matrix
J = edRef*dNr';
J1 = edN1*dNr';
%% Material parameters
a = obj.MAT.c1;
b = obj.MAT.c2;
c1 = obj.MAT.c;
d1 = obj.MAT.d;
c2 = obj.MAT.cc;
d2 = obj.MAT.dd;
e0 = obj.MAT.e0;
er = obj.MAT.e1;
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
thetaR = obj.MAT.thetaR;
k0 = obj.MAT.k0;
rho0 = obj.MAT.rhoSource;
%% Initialization
RX = rData{1,1};
RP = rData{2,1};
RD = rData{3,1};
RT = rData{4,1};
KXX = kData{1,1}; KXP = kData{1,2}; KXD = kData{1,3}; KXT = kData{1,4};
KPX = kData{2,1}; KPP = kData{2,2}; KPD = kData{2,3}; KPT = kData{2,4};
KDX = kData{3,1}; KDP = kData{3,2}; KDD = kData{3,3}; KDT = kData{3,4};
KTX = kData{4,1}; KTP = kData{4,2}; KTD = kData{4,3}; KTT = kData{4,4};
if stress ~= 0 || indexD ~= 0 % stress computation
    RX = zeros(size(N,1),1);    
    KXX = zeros(size(N,1),size(N,1));
end
GlobalEnergy.HelmholtzN1 = 0;
GlobalEnergy.EntropyN1 = 0;
GlobalEnergy.InternalEnergyN1 = 0;
GlobalEnergy.DeltaInternalEnergy = 0;  
%% DOFs
numberOfMechanicalDOFs = numel(RX);
numberOfElectricalDOFs = numel(RP);
numberOfInternalDOFs = numel(RD);
numberOfTemperatureDOFs = numel(RT);
% 
numberOfInternalNodes = size(obj.SHAPEF_D.N,2);
%% Gauss loop
for k = 1:NGP
    index = DIM*k-(DIM-1):DIM*k;
    detJ = det(J(:,index)');
    detJ1 = det(J1(:,index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,index)')\dNr(index,:);
    % Temperature
    thetaN1e = N(k,:)*thetaN1';
    thetaNe = N(k,:)*thetaN';
    DotTheta = (thetaN1e-thetaNe)/DT;  
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
    CxN1 = FxN1'*FxN1;
    CxN = FxN'*FxN;
    DotCx = (CxN1-CxN)/DT;
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN = 0.5*wedge(CxN,CxN);
    cxN1 = det(CxN1);
    cxN = det(CxN);
    % Electrical quantities
    EN1 = -dNx*phiN1';
    EN = -dNx*phiN';
    DN1 = zeros(3,1);
    for m = 1:numberOfInternalNodes
        DN1 = DN1 + N_D(k,m)*DN1e(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    DN = zeros(3,1);
    for m = 1:numberOfInternalNodes
        DN = DN + N_D(k,m)*DNe(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    % B-matrix (current configuration)
    BN1 = zeros(6,DIM*numel(obj.EDOF{j}));
    BN1(1,1:3:end) = FxN1(1,1)*dNx(1,:);
    BN1(1,2:3:end) = FxN1(2,1)*dNx(1,:);
    BN1(1,3:3:end) = FxN1(3,1)*dNx(1,:);
    BN1(2,1:3:end) = FxN1(1,2)*dNx(2,:);
    BN1(2,2:3:end) = FxN1(2,2)*dNx(2,:);
    BN1(2,3:3:end) = FxN1(3,2)*dNx(2,:);
    BN1(3,1:3:end) = FxN1(1,3)*dNx(3,:);
    BN1(3,2:3:end) = FxN1(2,3)*dNx(3,:);
    BN1(3,3:3:end) = FxN1(3,3)*dNx(3,:);
    BN1(4,1:3:end) = FxN1(1,1)*dNx(2,:)+FxN1(1,2)*dNx(1,:);
    BN1(4,2:3:end) = FxN1(2,1)*dNx(2,:)+FxN1(2,2)*dNx(1,:);
    BN1(4,3:3:end) = FxN1(3,1)*dNx(2,:)+FxN1(3,2)*dNx(1,:);
    BN1(5,1:3:end) = FxN1(1,2)*dNx(3,:)+FxN1(1,3)*dNx(2,:);
    BN1(5,2:3:end) = FxN1(2,2)*dNx(3,:)+FxN1(2,3)*dNx(2,:);
    BN1(5,3:3:end) = FxN1(3,2)*dNx(3,:)+FxN1(3,3)*dNx(2,:);
    BN1(6,1:3:end) = FxN1(1,1)*dNx(3,:)+FxN1(1,3)*dNx(1,:);
    BN1(6,2:3:end) = FxN1(2,1)*dNx(3,:)+FxN1(2,3)*dNx(1,:);
    BN1(6,3:3:end) = FxN1(3,1)*dNx(3,:)+FxN1(3,3)*dNx(1,:);
    BN1 = 2*BN1;
    %% energy
    PsiTN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiMIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiMVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
%     
    PsiTN = kappa*(thetaNe - thetaR - thetaNe*log(thetaNe/thetaR));
    PsiTMN = DIM*beta*c2*(thetaNe - thetaR)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiMIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiMVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
%     
    s0N1 = kappa*log(thetaN1e/thetaR) - DIM*beta*c2*(cxN1 - 1);
    s0N = kappa*log(thetaNe/thetaR) - DIM*beta*c2*(cxN - 1);
    u0N1 = PsiN1 + DN1'*EN1 - thetaN1e*s0N1;
    u0N = PsiN + DN'*EN - thetaNe*s0N;
%     
    GlobalEnergy.HelmholtzN1 = GlobalEnergy.HelmholtzN1 + PsiN1*detJ*wp(k);
    GlobalEnergy.EntropyN1 = GlobalEnergy.EntropyN1 + s0N1*detJ*wp(k);
    GlobalEnergy.InternalEnergyN1 = GlobalEnergy.InternalEnergyN1 + u0N1*detJ*wp(k);
    GlobalEnergy.DeltaInternalEnergy = GlobalEnergy.DeltaInternalEnergy + (u0N1 - u0N)*detJ*wp(k);
    %% residual
    DPsiEMN1_CN1 = 1/(2*er*e0*cxN1^(1/2))*DN1*DN1';
    DPsiMIsoN1_CN1 = a*eye(3);
    DPsi_CN1 = DPsiEMN1_CN1 + DPsiMIsoN1_CN1;
    DPsi_GN1 = b*eye(3);
    DPsiTMN1_cN1 = DIM*beta*c2*(thetaN1e - thetaR);
    DPsiEMN1_cN1 = -1/(4*er*e0*cxN1^(3/2))*DN1'*CxN1*DN1;
    DPsiMVolN1_cN1 = c1/2*(1-cxN1^(-1/2)) - d1/(2*cxN1);
    DPsi_cN1 = DPsiTMN1_cN1 + DPsiEMN1_cN1 + DPsiMVolN1_cN1;
    SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
    QxN1 = -k0*(1/cxN1*GxN1*(dNx*thetaN1'));
    DPsi_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;
    Du0_thetaN1 = kappa;
    Ds0_thetaN1 = kappa/thetaN1e;
    Ds0_cxN1 = -DIM*beta*c2;
%         
    if isa(rho0,'function_handle')
        XGP = kron(N(k,:),eye(DIM))*edRef(:);
        rho0GP = rho0(XGP(1),XGP(2),XGP(3));
    else
        rho0GP = rho0;
    end    
    if stress~=0 || indexD~=0 % von mises stress for plot
        if stress~=0
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            RX = RX + N(k,:)'*tempStress*detJ1*wp(k);
        elseif indexD~=0
            RX = RX + N(k,:)'*DN1(indexD)*detJ1*wp(k);
        end
        KXX = KXX + (N(k,:)'*N(k,:))*detJ*wp(k);
    else 
        %% Residual
        RX = RX + BN1'*0.5*SN1(map.Voigt)*detJ*wp(k);
        RP = RP + (dNx'*DN1 + obj.materialObjective.timeFunctionRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
        RD = RD + kron(N_D(k,:)',(DPsi_DN1 - EN1))*detJ*wp(k);        
        RT = RT + (N(k,:)'*DotTheta - (Du0_thetaN1)^(-1)*dNx'*QxN1 + N(k,:)'*Ds0_thetaN1^(-1)*Ds0_cxN1*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*map.Isym)')*detJ*wp(k);
%         RT = RT + (N(k,:)'*thetaN1e*(s0N1-s0N)/DT-dNx'*QxN1)*detJ*wp(k);
        %% Tangent
        % ---------------------------- KXX ----------------------------
        D2Psi_CN1cN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        D2Psi_cN1cN1 = c1/(4*cxN1^(3/2)) + d1/(2*cxN1^2) + 3/(8*er*e0)*cxN1^(-5/2)*DN1'*(CxN1*DN1);
        D2Psi_cN1CN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        elasticityTensor = 4*D2Psi_CN1cN1(map.Voigt)*GxN1(map.Voigt)' + 2*secondDiffFunction(DPsi_GN1) + 2*DPsi_cN1*secondDiffFunction(CxN1) + 4*D2Psi_cN1cN1*(GxN1(map.Voigt)*GxN1(map.Voigt)') + 4*GxN1(map.Voigt)*D2Psi_cN1CN1(map.Voigt)';
        materialTangent = 0.25*BN1'*elasticityTensor*BN1;
        % Geometrical part of tangent operator
        geometricalTangentPart = dNx'*SN1*dNx;
        geometricalTangent = zeros(numberOfMechanicalDOFs);
        for g = 1:DIM
            geometricalTangent(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = geometricalTangentPart;
        end
        KXX = KXX + (materialTangent + geometricalTangent)*detJ*wp(k);       
        % ---------------------------- KXP ----------------------------
        KXP = zeros(numberOfMechanicalDOFs,numberOfElectricalDOFs);
        % ---------------------------- KXD ----------------------------
        BEN1 = zeros(6,numberOfInternalDOFs);
        BEN1(1,1:3:end) = 2*DN1(1)*N_D(k,:);
        BEN1(2,2:3:end) = 2*DN1(2)*N_D(k,:);
        BEN1(3,3:3:end) = 2*DN1(3)*N_D(k,:);
        BEN1(4,1:3:end) = DN1(2)*N_D(k,:);
        BEN1(4,2:3:end) = DN1(1)*N_D(k,:);
        BEN1(5,2:3:end) = DN1(3)*N_D(k,:);
        BEN1(5,3:3:end) = DN1(2)*N_D(k,:);
        BEN1(6,1:3:end) = DN1(3)*N_D(k,:);
        BEN1(6,3:3:end) = DN1(1)*N_D(k,:);
        D2Psi_cN1DN1 = -1/(2*er*e0)*cxN1^(-3/2)*(CxN1*DN1);
        temp1 = BN1'*GxN1(map.Voigt);
        temp2 = kron(N_D(k,:),D2Psi_cN1DN1);
        KXD = KXD + (BN1'*0.5*2/(2*er*e0*cxN1^(1/2))*BEN1 + temp1*temp2(:)')*detJ*wp(k);        
        % ---------------------------- KXT ----------------------------
        D2Psi_cN1thetaN1 = DIM*beta*c2;
        DSN1_thetaN1 = 2*D2Psi_cN1thetaN1*GxN1;
        KXT = KXT + (BN1'*0.5*DSN1_thetaN1(map.Voigt))*N(k,:)*detJ*wp(k);        
        % ---------------------------- KPX ----------------------------
%       
        % ---------------------------- KPP ----------------------------
% 
        % ---------------------------- KPD ----------------------------
        KPD = KPD + kron(N_D(k,:),dNx')*detJ*wp(k);
        % ---------------------------- KPT ----------------------------
%         *Du0_thetaN1^(-1)
        % ---------------------------- KDX ----------------------------
        KDX = KDX + (BEN1'*0.5*2*1/(2*er*e0*cxN1^(1/2))*BN1 + temp2(:)*temp1')*detJ*wp(k);
        % ---------------------------- KDP ----------------------------
        KDP = KDP + kron(N_D(k,:),dNx')'*detJ*wp(k);
        % ---------------------------- KDD ----------------------------
        D2Psi_DN1DN1 = 1/(er*e0*cxN1^(1/2))*CxN1;
        KDD = KDD + kron(N_D(k,:)'*N_D(k,:),D2Psi_DN1DN1)*detJ*wp(k);     
        % ---------------------------- KDT ----------------------------
%         
        % ---------------------------- KTX ----------------------------
        D2s0_cxN1cxN1 = 0;
        kTX1 = N(k,:)'*(Ds0_thetaN1)^(-1)*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*map.Isym)'*D2s0_cxN1cxN1*(BN1'*GxN1(map.Voigt))';
        temp = wedge(CxN1,DotCx);
        kTX2 = N(k,:)'*(Ds0_thetaN1)^(-1)*Ds0_cxN1*(BN1'*temp(map.Voigt))';
        kTX3 = N(k,:)'*(Ds0_thetaN1)^(-1)*Ds0_cxN1*(1/DT*BN1'*GxN1(map.Voigt))';
        kTX4 = -dNx'*Du0_thetaN1^(-1)*k0*cxN1^(-2)*(GxN1*(dNx*thetaN1'))*(BN1'*GxN1(map.Voigt))';
        kTX5 = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for jj = 1:numberOfMechanicalDOFs
            BN1jj = map.Isyminv*BN1(:,jj);
            kTX5(:,jj) = dNx'*Du0_thetaN1^(-1)*k0*cxN1^(-1)*wedge(CxN1,BN1jj(map.VoigtInv))*(dNx*thetaN1');
        end
        KTX = KTX + (kTX1 + kTX2 + kTX3 + kTX4 + kTX5)*detJ*wp(k);
        % ---------------------------- KTP ----------------------------
%         
        % ---------------------------- KTD ----------------------------
%         
        % ---------------------------- KTT ----------------------------
        D2s0_thetaN1thetaN1 = -kappa*thetaN1e^(-2)*N(k,:);
        DQN1_thetaN1 = -k0*(1*cxN1^(-1)*GxN1*dNx*eye(size(dNx,2)));
        KTT = KTT + (N(k,:)'*1/DT*N(k,:) - Du0_thetaN1^(-1)*dNx'*DQN1_thetaN1 - N(k,:)'*Ds0_thetaN1^(-2)*D2s0_thetaN1thetaN1*(DotCx(map.Voigt)'*(GxN1(map.Voigt)'*map.Isym)'*Ds0_cxN1))*detJ*wp(k);
    end
end
if ~computePostData
    rData{1} = RX;
    rData{2} = RP;
    rData{3} = RT;
    rData{4} = RD;
    kData{1,1} = KXX;   kData{1,2} = KXP;   kData{1,3} = KXT;   kData{1,4} = KXD;
    kData{2,1} = KPX;   kData{2,2} = KPP;   kData{2,3} = KPT;   kData{2,4} = KPD;
    kData{3,1} = KTX;   kData{3,2} = KTP;   kData{3,3} = KTT;   kData{3,4} = KTD;
    kData{4,1} = KDX;   kData{4,2} = KDP;   kData{4,3} = KDT;   kData{4,4} = KDD;
end
end

function [rData,kData,GlobalEnergy] = Residuum3(j,DOFs,edofLocal,obj,map,time,DT,stress,indexD,rData,kData,dofsNT)
%% Shape functions
N = obj.SHAPEF.N;
N_D = obj.SHAPEF_D.N;
wp = obj.SHAPEF.wp;
NGP = obj.NGP;
dNr = obj.SHAPEF.dNr;
%% DOFs
DIM = obj.DIM;
% 
edN1 = dofsNT.edN1; 
phiN1 = dofsNT.phiN1;
DN1e = dofsNT.DN1;
thetaN1 = dofsNT.thetaN1;
% 
edN  = obj.QN(edofLocal,1:DIM)';
phiN = obj.QN(edofLocal,DIM+1)';
DNe  = obj.DN;
thetaN = obj.QN(edofLocal,DIM+2)';
% 
edRef = obj.QREF(edofLocal,1:DIM)';
%% Jacobi-Matrix
J = edRef*dNr';
J1 = edN1*dNr';
%% Material parameters
a = obj.MAT.c1;
b = obj.MAT.c2;
c1 = obj.MAT.c;
d1 = obj.MAT.d;
c2 = obj.MAT.cc;
d2 = obj.MAT.dd;
e0 = obj.MAT.e0;
er = obj.MAT.e1;
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
thetaR = obj.MAT.thetaR;
k0 = obj.MAT.k0;
rho0 = obj.MAT.rhoSource;
%% Initialization
RX = rData{1,1};
RP = rData{2,1};
RD = rData{3,1};
RT = rData{4,1};
KXX = kData{1,1}; KXP = kData{1,2}; KXD = kData{1,3}; KXT = kData{1,4};
KPX = kData{2,1}; KPP = kData{2,2}; KPD = kData{2,3}; KPT = kData{2,4};
KDX = kData{3,1}; KDP = kData{3,2}; KDD = kData{3,3}; KDT = kData{3,4};
KTX = kData{4,1}; KTP = kData{4,2}; KTD = kData{4,3}; KTT = kData{4,4};
if stress ~= 0 || indexD ~= 0 % stress computation
    RX = zeros(size(N,1),1);    
    KXX = zeros(size(N,1),size(N,1));
end
GlobalEnergy.HelmholtzN1 = 0;
GlobalEnergy.EntropyN1 = 0;
GlobalEnergy.InternalEnergyN1 = 0;
GlobalEnergy.DeltaInternalEnergy = 0;    
%% Gauss loop
for k = 1:NGP
    index = DIM*k-(DIM-1):DIM*k;
    detJ = det(J(:,index)');
    detJ1 = det(J1(:,index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,index)')\dNr(index,:);
    
    % Temperature
    thetaN1e = N(k,:)*thetaN1';
    thetaNe = N(k,:)*thetaN';
    DotTheta = (thetaN1e - thetaNe)/DT;
    
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
    CxN1 = FxN1'*FxN1;
    CxN = FxN'*FxN;
    DotCx = (CxN1-CxN)/DT;
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN = 0.5*wedge(CxN,CxN);
    DotGx = (GxN1-GxN)/DT;
    cxN1 = det(CxN1);
    cxN = det(CxN);
    Dotcx = (cxN1-cxN)/DT;
    
    % Electrical quantities
    EN1 = -dNx*phiN1';
    EN = -dNx*phiN';
    DN1 = zeros(3,1);
    for m = 1:DOFs.numNodes_D
        DN1 = DN1 + N_D(k,m)*DN1e(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    DN = zeros(3,1);
    for m = 1:DOFs.numNodes_D
        DN = DN + N_D(k,m)*DNe(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    
    % B-matrix (current configuration)
    BN1 = zeros(6,DOFs.numberOfMechanicalDOFs);
    BN1(1,1:3:end) = FxN1(1,1)*dNx(1,:);
    BN1(1,2:3:end) = FxN1(2,1)*dNx(1,:);
    BN1(1,3:3:end) = FxN1(3,1)*dNx(1,:);
    BN1(2,1:3:end) = FxN1(1,2)*dNx(2,:);
    BN1(2,2:3:end) = FxN1(2,2)*dNx(2,:);
    BN1(2,3:3:end) = FxN1(3,2)*dNx(2,:);
    BN1(3,1:3:end) = FxN1(1,3)*dNx(3,:);
    BN1(3,2:3:end) = FxN1(2,3)*dNx(3,:);
    BN1(3,3:3:end) = FxN1(3,3)*dNx(3,:);
    BN1(4,1:3:end) = FxN1(1,1)*dNx(2,:)+FxN1(1,2)*dNx(1,:);
    BN1(4,2:3:end) = FxN1(2,1)*dNx(2,:)+FxN1(2,2)*dNx(1,:);
    BN1(4,3:3:end) = FxN1(3,1)*dNx(2,:)+FxN1(3,2)*dNx(1,:);
    BN1(5,1:3:end) = FxN1(1,2)*dNx(3,:)+FxN1(1,3)*dNx(2,:);
    BN1(5,2:3:end) = FxN1(2,2)*dNx(3,:)+FxN1(2,3)*dNx(2,:);
    BN1(5,3:3:end) = FxN1(3,2)*dNx(3,:)+FxN1(3,3)*dNx(2,:);
    BN1(6,1:3:end) = FxN1(1,1)*dNx(3,:)+FxN1(1,3)*dNx(1,:);
    BN1(6,2:3:end) = FxN1(2,1)*dNx(3,:)+FxN1(2,3)*dNx(1,:);
    BN1(6,3:3:end) = FxN1(3,1)*dNx(3,:)+FxN1(3,3)*dNx(1,:);
    BN1 = 2*BN1;
    
    %% energy
    PsiThermoN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiThermoN1 + PsiTMN1 + thetaN1e/thetaR*(PsiEMN1 + PsiIsoN1 + PsiVolN1);
    
    PsiThermoN = kappa*(thetaNe - thetaR - thetaNe*log(thetaNe/thetaR));
    PsiTMN = DIM*beta*c2*(thetaNe - thetaR)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiThermoN + PsiTMN + thetaN1e/thetaR*(PsiEMN + PsiIsoN + PsiVolN);
    
    s0N1 = kappa*log(thetaN1e/thetaR) - DIM*beta*c2*(cxN1 - 1) - 1/thetaR*(PsiEMN1 + PsiIsoN1 + PsiVolN1);
    s0N = kappa*log(thetaNe/thetaR) - DIM*beta*c2*(cxN - 1) - 1/thetaR*(PsiEMN + PsiIsoN + PsiVolN);
    u0N1 = PsiN1 + DN1'*EN1 - thetaN1e*s0N1;
    u0N = PsiN + DN'*EN - thetaNe*s0N;
    
    GlobalEnergy.HelmholtzN1 = GlobalEnergy.HelmholtzN1 + PsiN1*detJ*wp(k);
    GlobalEnergy.EntropyN1 = GlobalEnergy.EntropyN1 + s0N1*detJ*wp(k);
    GlobalEnergy.InternalEnergyN1 = GlobalEnergy.InternalEnergyN1 + u0N1*detJ*wp(k);
    GlobalEnergy.DeltaInternalEnergy = GlobalEnergy.DeltaInternalEnergy + (u0N1 - u0N)*detJ*wp(k);
    
    %% residual
    DPsiEMN1_CN1 = 1/(2*er*e0*cxN1^(1/2))*DN1*DN1';
    DPsiIsoN1_CN1 = a*eye(3);
    DPsi_CN1 = thetaN1e/thetaNe*(DPsiEMN1_CN1 + DPsiIsoN1_CN1);
    DPsiIsoN1_GN1 = b*eye(3);
    DPsi_GN1 = thetaN1e/thetaNe*DPsiIsoN1_GN1;
    DPsiTMN1_cN1 = DIM*beta*(thetaN1e - thetaR)*(c2/2*cxN1^(-1/2) + d2*cxN1^(-3/2));
    DPsiEMN1_cN1 = -1/(4*er*e0*cxN1^(3/2))*DN1'*CxN1*DN1;
    DPsiVolN1_cN1 = c1/2*(1-cxN1^(-1/2)) - d1/(2*cxN1);
    DPsi_cN1 = DPsiTMN1_cN1 + thetaN1e/thetaNe*(DPsiEMN1_cN1 + DPsiVolN1_cN1);
    SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
    QxN1 = -k0*(1/cxN1*GxN1*(dNx*thetaN1'));
    DPsi_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;
    Du0_thetaN1 = kappa;
    Ds0_thetaN1 = kappa/thetaN1e;
    Ds0_cxN1 = -DIM*beta*c2;
    
    DPsiEMN1_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;

    Ds0_CxN1 = -1/thetaR*(DPsiEMN1_CN1 + DPsiIsoN1_CN1);
    Ds0_GxN1 = -1/thetaR*DPsiIsoN1_GN1;
    Ds0_DN1 = -1/thetaR*DPsiEMN1_DN1;
    
    if isa(rho0,'function_handle')
        XGP = kron(N(k,:),eye(DIM))*edREF(:);
        rho0GP = rho0(XGP(1),XGP(2),XGP(3));
    else
        rho0GP = rho0;
    end    
    if stress~=0 || indexD~=0 % von mises stress for plot
        if stress~=0
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            RX = RX + N(k,:)'*tempStress*detJ1*wp(k);
        elseif indexD~=0
            RX = RX + N(k,:)'*DN1(indexD)*detJ1*wp(k);
        end
        KXX = KXX + (N(k,:)'*N(k,:))*detJ*wp(k);
    else % residual and tangent
        RX = RX + BN1'*0.5*SN1(map.Voigt)*detJ*wp(k);
        RP = RP + (dNx'*DN1 + obj.MAT.timeFunctionRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
        RD = RD + kron(N_D(k,:)',(DPsi_DN1-EN1))*detJ*wp(k);
%         RT = RT + N(k,:)'*Du0_thetaN1/Ds0_thetaN1*Ds0_cxN1*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*Isym)' + N(k,:)'*Du0_thetaN1*DotTheta - dNx'*QxN1)'*detJ*wp(k);
        RT = RT + (N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT - (Du0_thetaN1)^(-1)*dNx'*QxN1 + N(k,:)'*Ds0_thetaN1^(-1)*Ds0_cxN1*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*map.Isym)')*detJ*wp(k);
%         RT = RT + (N(k,:)'*thetaN1e*(s0N1-s0N)/DT-dNx'*QxN1)*detJ*wp(k);       
%         RT = RT + (N(k,:)'*(Du0_thetaN1/Ds0_thetaN1*(Ds0_CxN1(map.Voigt)'*DotCx(map.Voigt) + Ds0_GxN1(map.Voigt)'*DotGx(map.Voigt) + Ds0_cxN1*Dotcx + Ds0_thetaN1*DotTheta + Ds0_DN1'*DN1)) - dNx'*QxN1)*detJ*wp(k);
    end
end
rData{1} = RX;
rData{2} = RP;
rData{3} = RD;
rData{4} = RT;
kData{1,1} = KXX;   kData{1,2} = KXP;   kData{1,3} = KXD;   kData{1,4} = KXT;
kData{2,1} = KPX;   kData{2,2} = KPP;   kData{2,3} = KPD;   kData{2,4} = KPT;
kData{3,1} = KDX;   kData{3,2} = KDP;   kData{3,3} = KDD;   kData{3,4} = KDT;
kData{4,1} = KTX;   kData{4,2} = KTP;   kData{4,3} = KTD;   kData{4,4} = KTT;
end

function [RX,RP,RD,RT,KXX,KXP,KXD,KXT,KPX,KPP,KPD,KPT,KDX,KDP,KDD,KDT,KTX,KTP,KTD,KTT,GlobalEnergy] = Residuum4(j,DOFs,edN,edN1,phiN,phiN1,DNe,DN1e,thetaN,thetaN1,QREF,QN1,edofLocal,dNr,objNT,map,time,DT,stress,indexD)
DIM = objNT.DIM;
edREF = QREF(edofLocal,1:DIM)';
% Run through all Gauss points
%% Shape functions
N = objNT.SHAPEF.N;
wp = objNT.SHAPEF.wp;
NGP = objNT.NGP;
%% Interpolation of mixed quantites
% Number of internal dofs
N_D = objNT.SHAPEF_D.N;
DOFs.numNodes_D = size(N_D,2);
DOFs.numberOfInternalDOFs = 3*DOFs.numNodes_D;
%% Material parameters
a = objNT.MAT.c1;
b = objNT.MAT.c2;
c1 = objNT.MAT.c;
d1 = objNT.MAT.d;
c2 = objNT.MAT.cc;
d2 = objNT.MAT.dd;
e0 = objNT.MAT.e0;
er = objNT.MAT.e1;
kappa = objNT.MAT.kappa;
beta = objNT.MAT.beta;
thetaR = objNT.MAT.thetaR;
k0 = objNT.MAT.k0;
rho0 = objNT.MAT.rhoSource;
    
% Initialize sub-residual vector
RX = zeros(DOFs.numberOfMechanicalDOFs,1);
RP = zeros(DOFs.numberOfElectricalDOFs,1);
RD = zeros(DOFs.numberOfInternalDOFs,1);
RT = zeros(DOFs.numberOfTemperatureDOFs,1);

% Initialize sub-tangent matrices
KXX = zeros(DOFs.numberOfMechanicalDOFs,DOFs.numberOfMechanicalDOFs);
KXP = zeros(DOFs.numberOfMechanicalDOFs,DOFs.numberOfElectricalDOFs);
KXD = zeros(DOFs.numberOfMechanicalDOFs,DOFs.numberOfInternalDOFs);
KXT = zeros(DOFs.numberOfMechanicalDOFs,DOFs.numberOfTemperatureDOFs);
KPX = zeros(DOFs.numberOfElectricalDOFs,DOFs.numberOfMechanicalDOFs);
KPP = zeros(DOFs.numberOfElectricalDOFs,DOFs.numberOfElectricalDOFs);
KPD = zeros(DOFs.numberOfElectricalDOFs,DOFs.numberOfInternalDOFs);
KPT = zeros(DOFs.numberOfElectricalDOFs,DOFs.numberOfTemperatureDOFs);
KDX = zeros(DOFs.numberOfInternalDOFs,DOFs.numberOfMechanicalDOFs);
KDP = zeros(DOFs.numberOfInternalDOFs,DOFs.numberOfElectricalDOFs);
KDD = zeros(DOFs.numberOfInternalDOFs,DOFs.numberOfInternalDOFs);
KDT = zeros(DOFs.numberOfInternalDOFs,DOFs.numberOfTemperatureDOFs);
KTX = zeros(DOFs.numberOfTemperatureDOFs,DOFs.numberOfMechanicalDOFs);
KTP = zeros(DOFs.numberOfTemperatureDOFs,DOFs.numberOfElectricalDOFs);
KTD = zeros(DOFs.numberOfTemperatureDOFs,DOFs.numberOfInternalDOFs);
KTT = zeros(DOFs.numberOfTemperatureDOFs,DOFs.numberOfTemperatureDOFs);

if stress ~= 0 || indexD ~= 0 % stress computation
    RX = zeros(size(N,1),1);    
    KXX = zeros(size(N,1),size(N,1));
end

GlobalEnergy.HelmholtzN1 = 0;
GlobalEnergy.EntropyN1 = 0;
GlobalEnergy.InternalEnergyN1 = 0;
GlobalEnergy.DeltaInternalEnergy = 0;
    
Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
Isyminv = [eye(3) zeros(3); zeros(3) 1/2*eye(3)];
J = QREF(edofLocal,1:DIM)'*dNr';
J1 = QN1(edofLocal,1:DIM)'*dNr';
% Run through all Gauss points
for k = 1:NGP
    indx = DIM*k-(DIM-1):DIM*k;
    detJ = det(J(:,indx)');
    detJ1 = det(J1(:,indx)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,indx)')\dNr(indx,:);
    
    % Temperature
    thetaN1e = N(k,:)*thetaN1';
    thetaNe = N(k,:)*thetaN';
    
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
    CxN1 = FxN1'*FxN1;
    CxN = FxN'*FxN;
    DotCx = 1/DT*(CxN1-CxN);
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN = 0.5*wedge(CxN,CxN);
    cxN1 = det(CxN1);
    cxN = det(CxN);
    
    % Electrical quantities
    EN1 = -dNx*phiN1';
    EN = -dNx*phiN';
    DN1 = zeros(3,1);
    for m = 1:DOFs.numNodes_D
        DN1 = DN1 + N_D(k,m)*DN1e(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    DN = zeros(3,1);
    for m = 1:DOFs.numNodes_D
        DN = DN + N_D(k,m)*DNe(j,ones(1,3)*(m - 1)*3 + (1:3))';
    end
    
    % B-matrix (current configuration)
    BN1 = zeros(6,DOFs.numberOfMechanicalDOFs);
    BN1(1,1:3:end) = FxN1(1,1)*dNx(1,:);
    BN1(1,2:3:end) = FxN1(2,1)*dNx(1,:);
    BN1(1,3:3:end) = FxN1(3,1)*dNx(1,:);
    BN1(2,1:3:end) = FxN1(1,2)*dNx(2,:);
    BN1(2,2:3:end) = FxN1(2,2)*dNx(2,:);
    BN1(2,3:3:end) = FxN1(3,2)*dNx(2,:);
    BN1(3,1:3:end) = FxN1(1,3)*dNx(3,:);
    BN1(3,2:3:end) = FxN1(2,3)*dNx(3,:);
    BN1(3,3:3:end) = FxN1(3,3)*dNx(3,:);
    BN1(4,1:3:end) = FxN1(1,1)*dNx(2,:)+FxN1(1,2)*dNx(1,:);
    BN1(4,2:3:end) = FxN1(2,1)*dNx(2,:)+FxN1(2,2)*dNx(1,:);
    BN1(4,3:3:end) = FxN1(3,1)*dNx(2,:)+FxN1(3,2)*dNx(1,:);
    BN1(5,1:3:end) = FxN1(1,2)*dNx(3,:)+FxN1(1,3)*dNx(2,:);
    BN1(5,2:3:end) = FxN1(2,2)*dNx(3,:)+FxN1(2,3)*dNx(2,:);
    BN1(5,3:3:end) = FxN1(3,2)*dNx(3,:)+FxN1(3,3)*dNx(2,:);
    BN1(6,1:3:end) = FxN1(1,1)*dNx(3,:)+FxN1(1,3)*dNx(1,:);
    BN1(6,2:3:end) = FxN1(2,1)*dNx(3,:)+FxN1(2,3)*dNx(1,:);
    BN1(6,3:3:end) = FxN1(3,1)*dNx(3,:)+FxN1(3,3)*dNx(1,:);
    BN1 = 2*BN1;
    
    %% energy
    PsiThermoN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiThermoN1 + PsiTMN1 + PsiEMN1 + PsiIsoN1 + PsiVolN1;
    
    PsiThermoN = kappa*(thetaNe - thetaR - thetaNe*log(thetaNe/thetaR));
    PsiTMN = DIM*beta*c2*(thetaNe - thetaR)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiThermoN + PsiTMN + PsiEMN + PsiIsoN + PsiVolN;
    
    s0N1 = kappa*log(thetaN1e/thetaR) - DIM*beta*c2*(cxN1 - 1);
    s0N = kappa*log(thetaNe/thetaR) - DIM*beta*c2*(cxN - 1);
    u0N1 = PsiN1 + DN1'*EN1 - thetaN1e*s0N1;
    u0N = PsiN + DN'*EN - thetaNe*s0N;
    
    GlobalEnergy.HelmholtzN1 = GlobalEnergy.HelmholtzN1 + PsiN1*detJ*wp(k);
    GlobalEnergy.EntropyN1 = GlobalEnergy.EntropyN1 + s0N1*detJ*wp(k);
    GlobalEnergy.InternalEnergyN1 = GlobalEnergy.InternalEnergyN1 + u0N1*detJ*wp(k);
    GlobalEnergy.DeltaInternalEnergy = GlobalEnergy.DeltaInternalEnergy + (u0N1 - u0N)*detJ*wp(k);
    
    %% residual
    DPsiIsoEMN1_CN1 = 1/(2*er*e0*cxN1^(1/2))*DN1*DN1';
    DPsiIsoN1_CN1 = a*eye(3);
    DPsi_CN1 = DPsiIsoEMN1_CN1 + DPsiIsoN1_CN1;
    DPsiIsoN1_GN1 = b*eye(3);
    DPsi_GN1 = DPsiIsoN1_GN1;
    DPsiTMN1_cN1 = DIM*beta*(thetaN1e - thetaR)*(c2/2*cxN1^(-1/2) + d2*cxN1^(-3/2));
    DPsiEMN1_cN1 = -1/(4*er*e0*cxN1^(3/2))*DN1'*CxN1*DN1;
    DPsiVolN1_cN1 = c1/2*(1-cxN1^(-1/2)) - d1/(2*cxN1);
    DPsi_cN1 = DPsiTMN1_cN1 + DPsiEMN1_cN1 + DPsiVolN1_cN1;
    SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
    QxN1 = -k0*(1/cxN1*GxN1*(dNx*thetaN1'));
    DPsi_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;
    Du0_thetaN1 = kappa;
    Ds0_thetaN1 = kappa/thetaN1e;
    Ds0_cxN1 = -DIM*beta*c2;
    
    if isa(rho0,'function_handle')
        XGP = kron(N(k,:),eye(DIM))*edREF(:);
        rho0GP = rho0(XGP(1),XGP(2),XGP(3));
    else
        rho0GP = rho0;
    end    
    if stress~=0 || indexD~=0 % von mises stress for plot
        if stress~=0
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            RX = RX + N(k,:)'*tempStress*detJ1*wp(k);
        elseif indexD~=0
            RX = RX + N(k,:)'*DN1(indexD)*detJ1*wp(k);
        end
        KXX = KXX + (N(k,:)'*N(k,:))*detJ*wp(k);
    else % residual and tangent
        RX = RX + BN1'*0.5*SN1(map.Voigt)*detJ*wp(k);
        RP = RP + (dNx'*DN1 + objNT.MAT.timeFktRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
        RD = RD + kron(N_D(k,:)',(DPsi_DN1-EN1))*detJ*wp(k);
        RT = RT + (N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT - (Du0_thetaN1)^(-1)*dNx'*QxN1 + N(k,:)'*Ds0_thetaN1^(-1)*Ds0_cxN1*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*Isym)')*detJ*wp(k);
%         RT = RT + (N(k,:)'*thetaN1e*(s0N1-s0N)/DT-dNx'*QxN1)*detJ*wp(k);        
    end
end
end

function secondDiffOperator = secondDiffFunction(D)
       secondDiffOperator = [   0          D(3,3)   D(2,2)      0           -D(3,2)     0;
                                D(3,3)     0        D(1,1)      0           0           -D(3,1);
                                D(2,2)     D(1,1)   0           -D(2,1)     0           0;
                                0          0        -D(2,1)     -D(3,3)     D(3,1)      D(2,3);
                                -D(3,2)	   0        0           D(3,1)      -D(1,1)     D(1,2);
                                0          -D(3,1)	0           D(3,2)      D(2,1)      -D(2,2) ];
        secondDiffOperator(1:3,1:6) = 2*secondDiffOperator(1:3,1:6);
        secondDiffOperator(4:6,1:3) = 2*secondDiffOperator(4:6,1:3);
end