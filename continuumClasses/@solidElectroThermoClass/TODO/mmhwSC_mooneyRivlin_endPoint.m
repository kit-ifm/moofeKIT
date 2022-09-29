function out = mmhwSC_mooneyRivlin_endPoint(obj,varargin)
% 05.04.2019 M.F.

%% Check input
DT = 1;
stress = 0;
indexD = 0;
time = 0;
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'stress')
        stress = varargin{i+1};
    elseif strcmp(varargin{i},'time')
        time = varargin{i+1};
    elseif strcmp(varargin{i},'d1')
        indexD = 1;
    elseif strcmp(varargin{i},'d2')
        indexD = 2;
    elseif strcmp(varargin{i},'d3')
        indexD = 3;
    end
end
flagNumericalTangent = obj.numericalTangent;
edof = obj.EDOF;
globalEdof = obj.GLOBALEDOF;
numberOfElements = numel(edof);
QN1 = obj.QN1;
DIM = obj.DIM;
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];
map.Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
map.Isyminv = [eye(3) zeros(3); zeros(3) 1/2*eye(3)];

out(numberOfElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
parfor j = 1:numberOfElements
% for j = 1:numberOfElements % for numerical tangent
    localEdof = edof{j};
    numDOFs = numel(globalEdof{j});
    if stress~=0 || indexD~=0
        [edofH1,edofH2] = expandEdof(globalEdof{j}(1:5:end));
        out(j).edofE = globalEdof{j}(1:5:end)';
        numStressDOFs = numDOFs/5;
        Re = zeros(numStressDOFs,1);
        Ke = zeros(numStressDOFs);
        stressDOFs = 1:numStressDOFs;
    else
        [edofH1,edofH2] = expandEdof(globalEdof{j});
        out(j).edofE = globalEdof{j}';
        Re = zeros(numDOFs,1);
        Ke = zeros(numDOFs);
    end
    out(j).pI = edofH1';
    out(j).pJ = edofH2';
%
    mechanicalDOFs = 1:numDOFs;
    mechanicalDOFs([4:5:numDOFs,5:5:numDOFs]) = [];
    electricalDOFs = 4:DIM+2:numDOFs;
    temperatureDOFs = 5:DIM+2:numDOFs;
    numberOfMechanicalDOFs = numel(mechanicalDOFs);
    numberOfElectricalDOFs = numel(electricalDOFs);
    numberOfInternalDOFs = DIM*size(obj.SHAPEF_D.N,2);
    numberOfTemperatureDOFs = numel(temperatureDOFs);
    numberOfFields = 4;     % X P D T
    dofsData = cell(numberOfFields,1);
    dofsData{1,1} = numberOfMechanicalDOFs;
    dofsData{2,1} = numberOfElectricalDOFs;
    dofsData{3,1} = numberOfInternalDOFs;
    dofsData{4,1} = numberOfTemperatureDOFs;
    rDataInitial = cell(numberOfFields,1);
    kDataInitial = cell(numberOfFields,numberOfFields);
    for ii = 1:numberOfFields
        rDataInitial{ii,1} = zeros(dofsData{ii,1},1);
        for jj = 1:numberOfFields
            kDataInitial{ii,jj} = zeros(dofsData{ii,1},dofsData{jj,1});
        end
    end
    if stress~=0 || indexD ~= 0
        %% stress computation
        [rData,kData,GlobalEnergy] = Residuum(j,obj,map,time,DT,stress,indexD,rDataInitial,kDataInitial,[],false);
        Re(stressDOFs,1) = rData{1,1};
        Ke(stressDOFs,stressDOFs) = kData{1,1};
    else
        %% residual and tangent computation
        [rData,kData,GlobalEnergy] = Residuum(j,obj,map,time,DT,stress,indexD,rDataInitial,kDataInitial,[],false);
%         % uncomment for numerical tangent
%         if flagNumericalTangent(1)
%             %% numerical tangent
%             % Deformation, electric potential and temperature
%             dofsNT.edN1 = QN1(localEdof,1:DIM)';
%             dofsNT.phiN1 = QN1(localEdof,DIM+1)';
%             dofsNT.DN1 = obj.DN1;
%             dofsNT.thetaN1 = QN1(localEdof,DIM+2)';
%             kDataNT = kDataInitial;
%             epsilonData = cell(numberOfFields,1);
%             for ii = 1:numberOfFields
%                 epsilonData{ii,1} = 1e-6;
%             end
%             for ii = 1:numberOfFields
%                 for jj = 1:numberOfFields
%                     for xx = 1:dofsData{jj,1}
%                         dofsNTxx = dofsNT;
%                         switch jj
%                             case 1  % X
%                                 dofsNTxx.edN1(xx) = dofsNTxx.edN1(xx) + epsilonData{jj,1};
%                             case 2  % P
%                                 dofsNTxx.phiN1(xx) = dofsNTxx.phiN1(xx) + epsilonData{jj,1};
%                             case 3  % D
%                                 dofsNTxx.DN1(j,xx) = dofsNTxx.DN1(j,xx) + epsilonData{jj,1};
%                             case 4  % T
%                                 dofsNTxx.thetaN1(xx) = dofsNTxx.thetaN1(xx) + epsilonData{jj,1};
%                             otherwise
%                                 error('number of DOFs exceeds proposed number');
%                         end
%                         [rxxData,~,~] = Residuum(j,obj,map,time,DT,stress,indexD,rDataInitial,kDataInitial,dofsNTxx,flagNumericalTangent(1));
%                         kDataNT{ii,jj}(:,xx) = (rxxData{ii,1}-rData{ii,1})/epsilonData{jj,1};
%                     end
%                     if flagNumericalTangent(2)
%                         disp(strcat('i:',num2str(ii),', j:', num2str(jj),'::',num2str(max(max(abs(kData{ii,jj} - kDataNT{ii,jj}))))))
%                     end
%                 end
%             end
%             kData = kDataNT;
%         end
        RX = rData{1};
        RP = rData{2};
        RD = rData{3};
        RT = rData{4};
        KXX = kData{1,1}; KXP = kData{1,2}; KXD = kData{1,3}; KXT = kData{1,4};
        KPX = kData{2,1}; KPP = kData{2,2}; KPD = kData{2,3}; KPT = kData{2,4};
        KDX = kData{3,1}; KDP = kData{3,2}; KDD = kData{3,3}; KDT = kData{3,4};
        KTX = kData{4,1}; KTP = kData{4,2}; KTD = kData{4,3}; KTT = kData{4,4};
%
        Re(mechanicalDOFs,1)= RX - KXD*(KDD\RD);
        Re(electricalDOFs,1)= RP - KPD*(KDD\RD);
        Re(temperatureDOFs,1)= RT - KTD*(KDD\RD);
%
        Ke(mechanicalDOFs,mechanicalDOFs) = KXX - KXD*(KDD\KDX);
        Ke(mechanicalDOFs,electricalDOFs) = KXP - KXD*(KDD\KDP);
        Ke(mechanicalDOFs,temperatureDOFs) = KXT - KXD*(KDD\KDT);
        Ke(electricalDOFs,mechanicalDOFs) = KPX - KPD*(KDD\KDX);
        Ke(electricalDOFs,electricalDOFs) = KPP - KPD*(KDD\KDP);
        Ke(electricalDOFs,temperatureDOFs) = KPT - KPD*(KDD\KDT);
        Ke(temperatureDOFs,mechanicalDOFs) = KTX - KTD*(KDD\KDX);
        Ke(temperatureDOFs,electricalDOFs) = KTP - KTD*(KDD\KDP);
        Ke(temperatureDOFs,temperatureDOFs) = KTT - KTD*(KDD\KDT);
%
        %% Save values and matrices (for update)
        zw(j).RD = rData{3,1};
        zw(j).KXD = kData{1,3};
        zw(j).KPD = kData{2,3};
        zw(j).KTD = kData{4,3};
        zw(j).KDX = kData{3,1}; zw(j).KDP = kData{3,2}; zw(j).KDT = kData{3,4}; zw(j).KDD = kData{3,3};
    end
    out(j).Re = Re;
    out(j).pK = Ke(:);
    out(j).ePot = [GlobalEnergy.InternalEnergy,...
                   GlobalEnergy.Helmholtz,...
                   GlobalEnergy.DE,...
                   GlobalEnergy.TS,...
                   GlobalEnergy.S,...
                   GlobalEnergy.DeltaS,...
                   GlobalEnergy.DeltaU];
end
if stress == 0 && indexD == 0
    obj.StaConStruct = zw;
end
end

function [rData,kData,GlobalEnergy] = Residuum(j,obj,map,time,DT,stress,indexD,rData,kData,dofsNT,flagNumericalTangent)
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
theta0 = obj.MAT.theta0;
k0 = obj.MAT.k0;
rho0 = obj.MAT.rhoSource;
R = obj.MAT.RSource;
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
    RX = zeros(size(N,2),1);
    KXX = zeros(size(N,2),size(N,2));
end
%
InternalEnergy = 0;
Helmholtz = 0;
DE = 0;
TS = 0;
S = 0;
DeltaS = 0;
DeltaU = 0;

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
    PsiTN = kappa*(thetaNe - theta0 - thetaNe*log(thetaNe/theta0));
    PsiTMN = -DIM*beta*c2*(thetaNe - theta0)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiMIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiMVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
    s0N = mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaNe,3);
    u0N = PsiN - DN'*EN + thetaNe*s0N;
    PsiTN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
%     PsiTMN1 = DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiTMN1 = -DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiMIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiMVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
    s0N1 = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,3);
    u0N1 = PsiN1 - DN1'*EN1 + thetaN1e*s0N1;
    % MR-Mike
%     Psi = a*(trace(CN1)-DIM) + b*(1/2*((trace(CN1))^2 - trace((CN1)^2))-DIM) - d1*log(sqrt(cxN1)) + kappa*(thetaN1 - theta0 -thetaN1*log(thetaN1/theta0)) + DIM*beta*(thetaN1 - theta0)*d1*I3N1^(-1);
    % MR-ThermoMech
%     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0)) - DIM*beta*(thetaN1e - theta0)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    % MR-ThermoElectroMech
%     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0)) - DIM*beta*(thetaN1e - theta0)*c2*(cxN1 - 1) + 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
%
    InternalEnergy = InternalEnergy + u0N1*detJ*wp(k);
    Helmholtz = Helmholtz + PsiN1*detJ*wp(k);
    DE = DE + DN1'*EN1*detJ*wp(k);
    TS = TS + thetaN1e*s0N1*detJ*wp(k);
    S = S + s0N1*detJ*wp(k);
    DeltaS = DeltaS + (s0N1-s0N)*detJ*wp(k);
    DeltaU = DeltaU + (u0N1-u0N)*detJ*wp(k);
    %% residual
    DPsi_CN1 = 1/(2*er*e0*cxN1^(1/2))*DN1*DN1' + a*eye(3);
    DPsi_GN1 = b*eye(3);
    DPsi_cN1 = -DIM*beta*c2*(thetaN1e - theta0) - 1/(4*er*e0*cxN1^(3/2))*DN1'*CxN1*DN1 + c1/2*(1-cxN1^(-1/2)) - d1/(2*cxN1);
    SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
    QxN1 = -k0*1/cxN1*GxN1*(dNx*thetaN1');
    DPsi_DN1 = 1/(er*e0*cxN1^(1/2))*CxN1*DN1;
    if isa(rho0,'function_handle')
        XGP = kron(N(k,:),eye(DIM))*edRef(:);
        rho0GP = rho0(XGP(1),XGP(2),XGP(3));
    else
        rho0GP = rho0;
    end
    if isa(R,'function_handle')
        XGP = kron(N(k,:),eye(DIM))*edRef(:);
        % FIXME!!!
%         RGP = R(XGP(1),XGP(2),XGP(3));
%           RGP = 0;
        RGP = R(XGP(2));
%         RGP = R(XGP(2));
        flagTransient = false;
    else
        flagTransient = true;
        RGP = R;
    end
    if ~(stress~=0 || indexD~=0) % von mises stress for plot
        %% Residual
        RX = RX + BN1'*0.5*SN1(map.Voigt)*detJ*wp(k);
        RP = RP + (dNx'*DN1 + obj.MAT.timeFktRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
        RD = RD + kron(N_D(k,:)',(DPsi_DN1 - EN1))*detJ*wp(k);
        if flagTransient
            RT = RT + (N(k,:)'*(thetaN1e*s0N1-thetaNe*s0N)/DT - N(k,:)'*(thetaN1e - thetaNe)/DT*s0N1 - dNx'*QxN1 - N(k,:)'*RGP)*detJ*wp(k);
%             RT = RT + (N(k,:)'*thetaNe*(s0N1 - s0N)/DT - dNx'*QxN1 + N(k,:)'*RGP)*detJ*wp(k);
        else
            RT = RT + (-dNx'*QxN1 - obj.MAT.timeFktRSource(time)*N(k,:)'*RGP)*detJ*wp(k);
        end
%         RT = RT + (N(k,:)'*thetaNe*(s0N1 - s0N)/DT - dNx'*QxN1)*detJ*wp(k);

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
        D2Psi_cN1thetaN1 = -DIM*beta*c2;
        DSN1_thetaN1 = 2*D2Psi_cN1thetaN1*GxN1;
        KXT = KXT + (BN1'*0.5*DSN1_thetaN1(map.Voigt))*N(k,:)*detJ*wp(k);
        % ---------------------------- KPX ----------------------------
%
        % ---------------------------- KPP ----------------------------
%
        % ---------------------------- KPD ----------------------------
        KPD = KPD + kron(N_D(k,:),dNx')*detJ*wp(k);
        % ---------------------------- KPT ----------------------------
%
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
        Ds0N1_cxN1 = DIM*beta*c2;
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
        KTX = KTX + (kTX1 + kTX2 + kTX3 + kTX4)*detJ*wp(k);
        % ---------------------------- KTP ----------------------------
%
        % ---------------------------- KTD ----------------------------
%
        % ---------------------------- KTT ----------------------------
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
        KTT = KTT + (kTT1 + kTT2 + kTT3)*detJ*wp(k);
    else
        if stress~=0
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            RX = RX + N(k,:)'*tempStress*detJ1*wp(k);
        elseif indexD~=0
            RX = RX + N(k,:)'*DN1(indexD)*detJ1*wp(k);
        end
        KXX = KXX + (N(k,:)'*N(k,:))*detJ*wp(k);
    end
end

GlobalEnergy.InternalEnergy = InternalEnergy;
GlobalEnergy.Helmholtz = Helmholtz;
GlobalEnergy.DE = DE;
GlobalEnergy.TS = TS;
GlobalEnergy.S = S;
GlobalEnergy.DeltaS = DeltaS;
GlobalEnergy.DeltaU = DeltaU;

rData{1} = RX;
rData{2} = RP;
rData{3} = RD;
rData{4} = RT;
kData{1,1} = KXX;   kData{1,2} = KXP;   kData{1,3} = KXD;   kData{1,4} = KXT;
kData{2,1} = KPX;   kData{2,2} = KPP;   kData{2,3} = KPD;   kData{2,4} = KPT;
kData{3,1} = KDX;   kData{3,2} = KDP;   kData{3,3} = KDD;   kData{3,4} = KDT;
kData{4,1} = KTX;   kData{4,2} = KTP;   kData{4,3} = KTD;   kData{4,4} = KTT;
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
theta0 = obj.MAT.theta0;
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
    PsiTN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiMIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiMVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
%
    PsiTN = kappa*(thetaNe - theta0 - thetaNe*log(thetaNe/theta0));
    PsiTMN = DIM*beta*c2*(thetaNe - theta0)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiMIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiMVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
%
    s0N1 = kappa*log(thetaN1e/theta0) - DIM*beta*c2*(cxN1 - 1);
    s0N = kappa*log(thetaNe/theta0) - DIM*beta*c2*(cxN - 1);
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
    DPsiTMN1_cN1 = DIM*beta*c2*(thetaN1e - theta0);
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
        RP = RP + (dNx'*DN1 + obj.MAT.timeFktRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
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
rData{1} = RX;
rData{2} = RP;
rData{3} = RD;
rData{4} = RT;
kData{1,1} = KXX;   kData{1,2} = KXP;   kData{1,3} = KXD;   kData{1,4} = KXT;
kData{2,1} = KPX;   kData{2,2} = KPP;   kData{2,3} = KPD;   kData{2,4} = KPT;
kData{3,1} = KDX;   kData{3,2} = KDP;   kData{3,3} = KDD;   kData{3,4} = KDT;
kData{4,1} = KTX;   kData{4,2} = KTP;   kData{4,3} = KTD;   kData{4,4} = KTT;
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
theta0 = obj.MAT.theta0;
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
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiThermoN1 + PsiTMN1 + thetaN1e/theta0*(PsiEMN1 + PsiIsoN1 + PsiVolN1);

    PsiThermoN = kappa*(thetaNe - theta0 - thetaNe*log(thetaNe/theta0));
    PsiTMN = DIM*beta*c2*(thetaNe - theta0)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiThermoN + PsiTMN + thetaN1e/theta0*(PsiEMN + PsiIsoN + PsiVolN);

    s0N1 = kappa*log(thetaN1e/theta0) - DIM*beta*c2*(cxN1 - 1) - 1/theta0*(PsiEMN1 + PsiIsoN1 + PsiVolN1);
    s0N = kappa*log(thetaNe/theta0) - DIM*beta*c2*(cxN - 1) - 1/theta0*(PsiEMN + PsiIsoN + PsiVolN);
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
    DPsiTMN1_cN1 = DIM*beta*(thetaN1e - theta0)*(c2/2*cxN1^(-1/2) + d2*cxN1^(-3/2));
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

    Ds0_CxN1 = -1/theta0*(DPsiEMN1_CN1 + DPsiIsoN1_CN1);
    Ds0_GxN1 = -1/theta0*DPsiIsoN1_GN1;
    Ds0_DN1 = -1/theta0*DPsiEMN1_DN1;

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
        RP = RP + (dNx'*DN1 + obj.MAT.timeFktRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
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
theta0 = objNT.MAT.theta0;
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
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiTMN1 = DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiThermoN1 + PsiTMN1 + PsiEMN1 + PsiIsoN1 + PsiVolN1;

    PsiThermoN = kappa*(thetaNe - theta0 - thetaNe*log(thetaNe/theta0));
    PsiTMN = DIM*beta*c2*(thetaNe - theta0)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiThermoN + PsiTMN + PsiEMN + PsiIsoN + PsiVolN;

    s0N1 = kappa*log(thetaN1e/theta0) - DIM*beta*c2*(cxN1 - 1);
    s0N = kappa*log(thetaNe/theta0) - DIM*beta*c2*(cxN - 1);
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
    DPsiTMN1_cN1 = DIM*beta*(thetaN1e - theta0)*(c2/2*cxN1^(-1/2) + d2*cxN1^(-3/2));
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
