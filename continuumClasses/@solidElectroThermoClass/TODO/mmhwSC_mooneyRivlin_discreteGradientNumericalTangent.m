function out = mmhwSC_mooneyRivlin_discreteGradientNumericalTangent(obj,varargin)
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
for j = 1:numberOfElements
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
    if stress~=0 || indexD ~= 0 
        %% stress computation
        [rData,kData,GlobalEnergy] = Residuum(j,obj,map,time,DT,stress,indexD,[],[],[],false);
        Re(stressDOFs,1) = rData{1,1};
        Ke(stressDOFs,stressDOFs) = kData{1,1};
    else
        %% residual and tangent computation
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
        flagDiscreteGradient = false(obj.NGP,4);
        [rData,kData,GlobalEnergy,flagDiscreteGradient] = Residuum(j,obj,map,time,DT,stress,indexD,rDataInitial,kDataInitial,[],false,flagDiscreteGradient); 
        if flagNumericalTangent(1)
            %% numerical tangent
            % Deformation, electric potential and temperature
            dofsNT.edN1 = QN1(localEdof,1:DIM)';
            dofsNT.phiN1 = QN1(localEdof,DIM+1)';
            dofsNT.DN1 = obj.DN1;
            dofsNT.thetaN1 = QN1(localEdof,DIM+2)';
            kDataNT = kDataInitial;
            epsilonData = cell(numberOfFields,1);
            for ii = 1:numberOfFields
                epsilonData{ii,1} = 1e-6;
            end
            for ii = 1:numberOfFields
                for jj = 1:numberOfFields
                    for xx = 1:dofsData{jj,1}
                        dofsNTxx = dofsNT;
                        switch jj
                            case 1  % X
                                dofsNTxx.edN1(xx) = dofsNTxx.edN1(xx) + epsilonData{jj,1};
                            case 2  % P
                                dofsNTxx.phiN1(xx) = dofsNTxx.phiN1(xx) + epsilonData{jj,1};
                            case 3  % D
                                dofsNTxx.DN1(j,xx) = dofsNTxx.DN1(j,xx) + epsilonData{jj,1};
                            case 4  % T
                                dofsNTxx.thetaN1(xx) = dofsNTxx.thetaN1(xx) + epsilonData{jj,1};
                            otherwise
                                error('number of DOFs exceeds proposed number');
                        end
                        [rxxData,~,~,~] = Residuum(j,obj,map,time,DT,stress,indexD,rDataInitial,kDataInitial,dofsNTxx,flagNumericalTangent(1),flagDiscreteGradient);
                        kDataNT{ii,jj}(:,xx) = (rxxData{ii,1}-rData{ii,1})/epsilonData{jj,1};
                    end
                    if flagNumericalTangent(2)
                        disp(strcat('i:',num2str(ii),', j:', num2str(jj),'::',num2str(max(max(abs(kData{ii,jj} - kDataNT{ii,jj}))))))
                    end
                end
            end
            kData = kDataNT;
        end
%       
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

function [rData,kData,GlobalEnergy,flagDiscreteGradient] = Residuum(j,obj,map,time,DT,stress,indexD,rData,kData,dofsNT,flagNumericalTangent,flagDiscreteGradient)
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
edN05 = 0.5*(edN + edN1);
phiN = obj.QN(edofLocal,DIM+1)';
phiN05 = 0.5*(phiN + phiN1);
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
    % Quantities at N05
    FxN05 = edN05*dNx';
    CxN05 = FxN05'*FxN05;
    GxN05 = 0.5*wedge(CxN05,CxN05);
    cxN05 = det(CxN05);
    EN05 = -dNx*phiN05';
    DN05 = 0.5*(DN + DN1);
    thetaN05e = 0.5*(thetaN1e + thetaNe);
    thetaN05 = 0.5*(thetaN1 + thetaN);
    % B-matrix (midpoint configuration)
    BN05 = zeros(6,DIM*numel(obj.EDOF{j}));
    BN05(1,1:3:end) = FxN05(1,1)*dNx(1,:);
    BN05(1,2:3:end) = FxN05(2,1)*dNx(1,:);
    BN05(1,3:3:end) = FxN05(3,1)*dNx(1,:);
    BN05(2,1:3:end) = FxN05(1,2)*dNx(2,:);
    BN05(2,2:3:end) = FxN05(2,2)*dNx(2,:);
    BN05(2,3:3:end) = FxN05(3,2)*dNx(2,:);
    BN05(3,1:3:end) = FxN05(1,3)*dNx(3,:);
    BN05(3,2:3:end) = FxN05(2,3)*dNx(3,:);
    BN05(3,3:3:end) = FxN05(3,3)*dNx(3,:);
    BN05(4,1:3:end) = FxN05(1,1)*dNx(2,:) + FxN05(1,2)*dNx(1,:);
    BN05(4,2:3:end) = FxN05(2,1)*dNx(2,:) + FxN05(2,2)*dNx(1,:);
    BN05(4,3:3:end) = FxN05(3,1)*dNx(2,:) + FxN05(3,2)*dNx(1,:);
    BN05(5,1:3:end) = FxN05(1,2)*dNx(3,:) + FxN05(1,3)*dNx(2,:);
    BN05(5,2:3:end) = FxN05(2,2)*dNx(3,:) + FxN05(2,3)*dNx(2,:);
    BN05(5,3:3:end) = FxN05(3,2)*dNx(3,:) + FxN05(3,3)*dNx(2,:);
    BN05(6,1:3:end) = FxN05(1,1)*dNx(3,:) + FxN05(1,3)*dNx(1,:);
    BN05(6,2:3:end) = FxN05(2,1)*dNx(3,:) + FxN05(2,3)*dNx(1,:);
    BN05(6,3:3:end) = FxN05(3,1)*dNx(3,:) + FxN05(3,3)*dNx(1,:);
    BN05 = 2*BN05;
    
    %% residual
    QxN05 = -k0*(1/cxN05*GxN05*(dNx*thetaN05'));
    s0N1 = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,3);
    s0N = mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaNe,3);
    % algorithmic evaluations
    GxAlgo = 0.5*(GxN1 + GxN);
    CxAlgo = 0.5*(CxN + CxN1);    
    GxAlgo2 = 1/3*(wedge(CxAlgo,CxAlgo) + GxAlgo);
    %% energy  
    PsiTN = kappa*(thetaNe - theta0 - thetaNe*log(thetaNe/theta0));
    PsiTMN = -DIM*beta*c2*(thetaNe - theta0)*(cxN - 1);
    PsiEMN = 1/(2*er*e0*cxN^(1/2))*DN'*CxN*DN;
    PsiMIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiMVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
    u0N = PsiN - DN'*EN + thetaNe*s0N;
    PsiTN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiTMN1 = -DIM*beta*c2*(thetaN1e - theta0)*(cxN1 - 1);
    PsiEMN1 = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    PsiMIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiMVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;    
    u0N1 = PsiN1 - DN1'*EN1 + thetaN1e*s0N1;
%     
    InternalEnergy = InternalEnergy + u0N1*detJ*wp(k);
    Helmholtz = Helmholtz + PsiN1*detJ*wp(k);
    DE = DE + DN1'*EN1*detJ*wp(k);
    TS = TS + thetaN1e*s0N1*detJ*wp(k);
    S = S + s0N1*detJ*wp(k);
    DeltaS = DeltaS + (s0N1-s0N)*detJ*wp(k);
    DeltaU = DeltaU + (u0N1-u0N)*detJ*wp(k);
    
    deltaC = CxN1-CxN;
    if (norm(deltaC) >= 10^(-8) && ~flagNumericalTangent) || (flagNumericalTangent && flagDiscreteGradient(k,1))
        if ~flagNumericalTangent
            flagDiscreteGradient(k,1) = true;
        end
        DPsi__C_CxN05_GN_cN_DN_TN = mooneyRivlin(obj,CxN05,GxN,cxN,DN,thetaNe,14);
        DPsi__C_CxN05_GN1_cN1_DN1_TN1 = mooneyRivlin(obj,CxN05,GxN1,cxN1,DN1,thetaN1e,14);
        DPsi__C = 0.5*(  DPsi__C_CxN05_GN_cN_DN_TN + (mooneyRivlin(obj,CxN1,GxN,cxN,DN,thetaNe,1) - mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaNe,1) - DPsi__C_CxN05_GN_cN_DN_TN(:)'*deltaC(:))/(deltaC(:)'*deltaC(:))*deltaC...
            + DPsi__C_CxN05_GN1_cN1_DN1_TN1 + (mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,1) - mooneyRivlin(obj,CxN,GxN1,cxN1,DN1,thetaN1e,1) - DPsi__C_CxN05_GN1_cN1_DN1_TN1(:)'*deltaC(:))/(deltaC(:)'*deltaC(:))*deltaC);
    else
        DPsi__C = mooneyRivlin(obj,CxN05,GxN05,cxN05,DN05,thetaN05e,14);
    end
    %       DPsi_G
    DPsi_GxN05 = b*eye(3);
    DPsi__G = DPsi_GxN05;
    %       DPsi_c
    deltac = cxN1-cxN;
    if (norm(deltac) >= 10^(-8) && ~flagNumericalTangent) || (flagNumericalTangent && flagDiscreteGradient(k,2))
        if ~flagNumericalTangent
            flagDiscreteGradient(k,2) = true;
        end
        DPsi__c = 0.5*(  (mooneyRivlin(obj,CxN1,GxN1,cxN1,DN,thetaNe,1) - mooneyRivlin(obj,CxN1,GxN1,cxN,DN,thetaNe,1))/deltac...
                      + (mooneyRivlin(obj,CxN,GxN,cxN1,DN1,thetaN1e,1) - mooneyRivlin(obj,CxN,GxN,cxN,DN1,thetaN1e,1))/deltac);
    else
        DPsi__c = mooneyRivlin(obj,CxN05,GxN05,cxN05,DN05,thetaN05e,15);
    end
    %       DPsi_D
    deltaD = DN1-DN;
    if (norm(deltaD) >= 10^(-8) && ~flagNumericalTangent) || (flagNumericalTangent && flagDiscreteGradient(k,3))
        if ~flagNumericalTangent
            flagDiscreteGradient(k,3) = true;
        end
        DPsi__D_CN1_GN1_cN1_DN05_TN = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN05,thetaNe,16);
        DPsi__D_CN_GN_cN_DN05_TN1 = mooneyRivlin(obj,CxN,GxN,cxN,DN05,thetaN1e,16);
        DPsi__D = 0.5*(  DPsi__D_CN1_GN1_cN1_DN05_TN + (mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaNe,1) - mooneyRivlin(obj,CxN1,GxN1,cxN1,DN,thetaNe,1) - DPsi__D_CN1_GN1_cN1_DN05_TN(:)'*deltaD(:))/(deltaD(:)'*deltaD(:))*deltaD...
            + DPsi__D_CN_GN_cN_DN05_TN1 + (mooneyRivlin(obj,CxN,GxN,cxN,DN1,thetaN1e,1) - mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaN1e,1) - DPsi__D_CN_GN_cN_DN05_TN1(:)'*deltaD(:))/(deltaD(:)'*deltaD(:))*deltaD);
    else
        DPsi__D = mooneyRivlin(obj,CxN05,GxN05,cxN05,DN05,thetaN05e,16);
    end
    %       DPsi_theta
    deltat = thetaN1e-thetaNe;
    if (norm(deltat) >= 10^(-8) && ~flagNumericalTangent) || (flagNumericalTangent && flagDiscreteGradient(k,4))
        if ~flagNumericalTangent
            flagDiscreteGradient(k,4) = true;
        end
        DPsi__t = 0.5*(  (mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,1) - mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaNe,1))/deltat...
                      + (mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaN1e,1) - mooneyRivlin(obj,CxN,GxN,cxN,DN,thetaNe,1))/deltat);
    else
        DPsi__t = mooneyRivlin(obj,CxN05,GxN05,cxN05,DN05,thetaN05e,17);
    end
    
    sAlgo = -DPsi__t;
    SAlgo = 2*(DPsi__C + wedge(DPsi__G,CxAlgo) + DPsi__c*GxAlgo2);
    if stress~=0 || indexD~=0 % von mises stress for plot
        if stress~=0
            DPsi_CN1 = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,14);
            DPsi_GN1 = b*eye(3);
            DPsi_cN1 = mooneyRivlin(obj,CxN1,GxN1,cxN1,DN1,thetaN1e,15);
            SN1 = 2*(DPsi_CN1 + wedge(DPsi_GN1,CxN1) + DPsi_cN1*GxN1);
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            RX = RX + N(k,:)'*tempStress*detJ1*wp(k);
        elseif indexD~=0
            RX = RX + N(k,:)'*DN1(indexD)*detJ1*wp(k);
        end
        KXX = KXX + (N(k,:)'*N(k,:))*detJ*wp(k);
    else 
        %% Residual
        RX = RX + BN05'*0.5*SAlgo(map.Voigt)*detJ*wp(k);
        RP = RP + dNx'*DN05*detJ*wp(k);         %  + obj.MAT.timeFktRhoSource(time)*N(k,:)'*rho0GP)*detJ*wp(k);
        RD = RD + kron(N_D(k,:)',(DPsi__D - EN05))*detJ*wp(k);
        RT = RT + (N(k,:)'*(thetaN1e*s0N1-thetaNe*s0N)/DT - N(k,:)'*(thetaN1e - thetaNe)/DT*sAlgo - dNx'*QxN05)*detJ*wp(k);

%         %% Tangent
%         % ---------------------------- KXX ----------------------------
%         D2Psi_CxN05CxN05 = -1/(4*er*e0)*cxN05^(-3/2)*DN05*DN05';
%         D2Psi_CxN05CxN05 = c1/(4*cxN05^(3/2)) + d1/(2*cxN05^2) + 3/(8*er*e0)*cxN05^(-5/2)*DN05'*(CxN05*DN05);
%         D2Psi_CxN05CxN05 = -1/(4*er*e0)*cxN05^(-3/2)*DN05*DN05';
%         elasticityTensor = 4*D2Psi_CxN05CxN05(map.Voigt)*GxN05(map.Voigt)' + 2*secondDiffFunction(DPsi_GxN05) + 2*DPsi_CxN05*secondDiffFunction(CxN05) + 4*D2Psi_CxN05CxN05*(GxN05(map.Voigt)*GxN05(map.Voigt)') + 4*GxN05(map.Voigt)*D2Psi_CxN05CxN05(map.Voigt)';
%         materialTangent = 0.25*BN05'*elasticityTensor*BN05;
%         % Geometrical part of tangent operator
%         geometricalTangentPart = dNx'*SN05*dNx;
%         geometricalTangent = zeros(numberOfMechanicalDOFs);
%         for g = 1:DIM
%             geometricalTangent(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = geometricalTangentPart;
%         end
%         KXX = KXX + 0.5*(materialTangent + geometricalTangent)*detJ*wp(k);       
%         % ---------------------------- KXP ----------------------------
%         KXP = zeros(numberOfMechanicalDOFs,numberOfElectricalDOFs);
%         % ---------------------------- KXD ----------------------------
%         BEN05 = zeros(6,numberOfInternalDOFs);
%         BEN05(1,1:3:end) = 2*DN05(1)*N_D(k,:);
%         BEN05(2,2:3:end) = 2*DN05(2)*N_D(k,:);
%         BEN05(3,3:3:end) = 2*DN05(3)*N_D(k,:);
%         BEN05(4,1:3:end) = DN05(2)*N_D(k,:);
%         BEN05(4,2:3:end) = DN05(1)*N_D(k,:);
%         BEN05(5,2:3:end) = DN05(3)*N_D(k,:);
%         BEN05(5,3:3:end) = DN05(2)*N_D(k,:);
%         BEN05(6,1:3:end) = DN05(3)*N_D(k,:);
%         BEN05(6,3:3:end) = DN05(1)*N_D(k,:); 
%         D2Psi_CxN05DN05 = -1/(2*er*e0)*cxN05^(-3/2)*(CxN05*DN05);
%         temp1 = BN05'*GxN05(map.Voigt);
%         temp2 = kron(N_D(k,:),D2Psi_CxN05DN05);
%         KXD = KXD + 0.5*(BN05'*0.5*2/(2*er*e0*cxN05^(1/2))*BEN05 + temp1*temp2(:)')*detJ*wp(k);        
%         % ---------------------------- KXT ----------------------------
%         D2Psi_CxN05thetaN05 = DIM*beta*c2;
%         DSN05_thetaN05 = 2*D2Psi_CxN05thetaN05*GxN05;
%         KXT = KXT + 0.5*(BN05'*DSN05_thetaN05(map.Voigt))*N(k,:)*detJ*wp(k);        
%         % ---------------------------- KPX ----------------------------
% %       
%         % ---------------------------- KPP ----------------------------
% % 
%         % ---------------------------- KPD ----------------------------
%         KPD = KPD + 0.5*kron(N_D(k,:),dNx')*detJ*wp(k);
%         % ---------------------------- KPT ----------------------------
% %         
%         % ---------------------------- KDX ----------------------------
%         KDX = KDX + 0.5*(BEN05'*0.5*2*1/(2*er*e0*cxN05^(1/2))*BN05 + temp2(:)*temp1')*detJ*wp(k);
%         % ---------------------------- KDP ----------------------------
%         KDP = KDP + 0.5*kron(N_D(k,:),dNx')'*detJ*wp(k);
%         % ---------------------------- KDD ----------------------------
%         D2Psi_DN05DN05 = 1/(er*e0*cxN05^(1/2))*CxN05;
%         KDD = KDD + 0.5*kron(N_D(k,:)'*N_D(k,:),D2Psi_DN05DN05)*detJ*wp(k);     
%         % ---------------------------- KDT ----------------------------
% %         
%         % ---------------------------- KTX ----------------------------
%         D2s0_cxN05cxN05 = 0;
%         kTX1 = N(k,:)'*(Ds0_thetaN05)^(-1)*DotCx(map.Voigt)'*(GxN05(map.Voigt)'*map.Isym)'*D2s0_cxN05cxN05*(BN05'*GxN05(map.Voigt))';
%         temp = wedge(CxN05,DotCx);
%         kTX2 = N(k,:)'*(Ds0_thetaN05)^(-1)*Ds0_cxN05*(BN05'*temp(map.Voigt))';
%         kTX3 = N(k,:)'*(Ds0_thetaN05)^(-1)*Ds0_cxN05*(1/DT*BN05'*GxN05(map.Voigt))';
%         kTX4 = -dNx'*Du0__theta^(-1)*k0*cxN05^(-2)*(GxN05*(dNx*thetaN05'))*(BN05'*GxN05(map.Voigt))';
%         kTX5 = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
%         for jj = 1:numberOfMechanicalDOFs
%             BN05jj = map.Isyminv*BN05(:,jj);
%             kTX5(:,jj) = dNx'*Du0__theta^(-1)*k0*cxN05^(-1)*wedge(CxN05,BN05jj(map.VoigtInv))*(dNx*thetaN05');
%         end
%         KTX = KTX + (0.5*kTX1 + 0.5*kTX2 + kTX3 + 0.5*kTX4 + 0.5*kTX5)*detJ*wp(k);
%         % ---------------------------- KTP ----------------------------
% %         
%         % ---------------------------- KTD ----------------------------
% %         
%         % ---------------------------- KTT ----------------------------
%         D2s0_thetaN05thetaN05 = -kappa*thetaN05e^(-2)*N(k,:);
%         DQN05_thetaN05 = -k0*(1*cxN05^(-1)*GxN05*dNx*eye(size(dNx,2)));
%         KTT = KTT + (N(k,:)'*1/DT*N(k,:) - 0.5*Du0__theta^(-1)*dNx'*DQN05_thetaN05 - 0.5*N(k,:)'*Ds0_thetaN05^(-2)*D2s0_thetaN05thetaN05*(DotCx(map.Voigt)'*(GxN05(map.Voigt)'*map.Isym)'*Ds0_cxN05))*detJ*wp(k);
    end
end
% 
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

function out = mooneyRivlin(obj,C,G,c,D,theta,choice)
a = obj.MAT.c1;
b = obj.MAT.c2;
c1 = obj.MAT.c;
d1 = obj.MAT.d;
c2 = obj.MAT.cc;
e0 = obj.MAT.e0;
er = obj.MAT.e1;
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
theta0 = obj.MAT.theta0;
DIM = obj.DIM;

PsiT = kappa*(theta - theta0 - theta*log(theta/theta0));
PsiTM = DIM*beta*c2*(theta - theta0)*(c - 1);
PsiEM = 1/(2*er*e0*c^(1/2))*D'*C*D;
PsiMIso = a*(trace(C)-3) + b*(trace(G)-3);
PsiMVol = c1/2*(sqrt(c)-1)^2 - d1*log(sqrt(c));
Psi = PsiT + PsiTM + PsiEM + PsiMIso + PsiMVol;
DPsi_C = 1/(2*er*e0*c^(1/2))*D*D' + a*eye(3);
DPsi_c = DIM*beta*c2*(theta - theta0) - 1/(4*er*e0*c^(3/2))*D'*C*D + c1/2*(1-c^(-1/2)) - d1/(2*c);
DPsi_D = 1/(er*e0*c^(1/2))*C*D;
DPsi_t = -kappa*log(theta/theta0) + DIM*beta*c2*(c - 1);


s0 = -DPsi_t;
Ds0_c = -DIM*beta*c2;

% u0 = Psi + theta*s0;
% u0 = Psi + theta*s0 - D'*E;
u0 = kappa*(theta-theta0) - DIM*beta*c2*theta0*(c-1) + 1/(2*er*e0*c^(1/2))*D'*C*D + a*(trace(C)-3) + b*trace(G) - d1*log(c^(1/2)) + c1/2*(c^(1/2) - 1)^2;
Du0_C = 1/(2*er*e0*c^(1/2))*D*D' + a*eye(3);
Du0_c = -DIM*beta*c2*theta0 - 1/(4*er*e0*c^(3/2))*D'*C*D + c1/2*(1 - c^(-1/2)) - d1/(2*c);
Du0_D = 1/(er*e0*c^(1/2))*C*D;

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
end
end

function secondDiffOperator = secondDiffFunction(D)
secondDiffOperator = [  0          D(3,3)   D(2,2)      0           -D(3,2)     0;
                        D(3,3)     0        D(1,1)      0           0           -D(3,1);
                        D(2,2)     D(1,1)   0           -D(2,1)     0           0;
                        0          0        -D(2,1)     -D(3,3)     D(3,1)      D(2,3);
                        -D(3,2)	   0        0           D(3,1)      -D(1,1)     D(1,2);
                        0          -D(3,1)	0           D(3,2)      D(2,1)      -D(2,2) ];
secondDiffOperator(1:3,1:6) = 2*secondDiffOperator(1:3,1:6);
secondDiffOperator(4:6,1:3) = 2*secondDiffOperator(4:6,1:3);
end
