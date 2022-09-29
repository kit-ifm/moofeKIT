function out = mmhwGenericCascadeSC_mooneyRivlin_endPoint(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = mmhwGenericCascadeSC_mooneyRivlinThermo_endPoint(obj,'PropertyName',PropertyValue)
%
% Description
% mmhwCascadeSC=: mixed method based with Hu-Washizu variational principle
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 30.05.2017 M.S., A.J., M.F.

%% Check input
DT = 1;
LocalDofs = [];
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'LocalDofs')
        LocalDofs = varargin{i+1};
    end
end
ansatzfunction = obj.Ansatzfunction;
dNrAll = obj.SHAPEF.dNr;
% Element stuff
%% Create residual and tangent
edof = obj.EDOF;
numElements = numel(edof);
edof = obj.GLOBALEDOF;
edofZw = obj.EDOF;
QREF = obj.QREF;
QN1 = obj.QN1;
QN = obj.QN;
% Internal variables N+1
SIGMA_CN1e  = obj.SIGMA_FN1;
SIGMA_GN1e  = obj.SIGMA_HN1;
SIGMA_cN1e  = obj.SIGMA_JN1;
CN1e  = obj.FN1;
GN1e  = obj.HN1;
cN1e  = obj.JN1e;
% Internal variables N
CNe  = obj.FN;
cNe  = obj.JNe;
% Element stuff
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];
Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
Isyminv = [eye(3) zeros(3); zeros(3) 1/2*eye(3)];
DIM = obj.DIM;
%% Shape functions
N = obj.SHAPEF.N;
wp = obj.SHAPEF.wp;
NGP = obj.NGP;
%% Material parameters
a = obj.MAT.c1;
b = obj.MAT.c2;
c = obj.MAT.c;
d = 2*(a+2*b);
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
theta0 = obj.MAT.theta0;
k0 = obj.MAT.k0;
%% Interpolation of mixed quantites
% Number of internal dofs
N_C = obj.SHAPEF_F.N;
N_G = obj.SHAPEF_Cof.N;
N_c =  obj.SHAPEF_Det.N;
numNodes_C = size(N_C,2);
numNodes_G = size(N_G,2);
numNodes_c = size(N_c,2);
numInternalDOFs = (numNodes_C+numNodes_G)*6+numNodes_c;
%% Element Output
out(numElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
zw  = struct('RSe',[],'KSe',[]);
parfor j = 1:numElements
    dNr = [];
    edofLocal = edofZw{j};
    [edofH1,edofH2] = expandEdof(edof{j});
    out(j).pI = edofH1';
    out(j).pJ = edofH2';
    out(j).edofE = edof{j}';
    numDOFs = numel(edof{j});
    if strcmpi(ansatzfunction,'Lagrange')
        dNr = dNrAll;
    elseif strcmpi(ansatzfunction,'HNURBS')
        dNr = dNrAll{j};
    end
    % Deformation and temperature
    edN1 = QN1(edofLocal,1:DIM)';
    edN = QN(edofLocal,1:DIM)';
    thetaN1 = QN1(edofLocal,DIM+1)';
    thetaN = QN(edofLocal,DIM+1)';
    % External dofs
    temperatureDOFs = 4:DIM+1:numDOFs;
    mechanicalDOFs = 1:numDOFs;
    mechanicalDOFs(4:DIM+1:numDOFs) = [];
    % Number of external dofs
    numberOfTemperatureDOFs = numel(temperatureDOFs);
    numberOfMechanicalDOFs = numel(mechanicalDOFs);
    %Initialize element resdiual and tangent matrix
    Re = zeros(numDOFs,1);
    Ke = zeros(numDOFs);
    %Initialize sub-tangent matrices
    KXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
    KXSIG = zeros(numberOfMechanicalDOFs,numInternalDOFs);
    KEPSEPS = zeros(numInternalDOFs,numInternalDOFs);
    KEPSSIG = zeros(numInternalDOFs,numInternalDOFs);
    KEPST = zeros(numInternalDOFs,numberOfTemperatureDOFs);
    KSIGX = zeros(numInternalDOFs,numberOfMechanicalDOFs);
    KSIGEPS = zeros(numInternalDOFs,numInternalDOFs);
    KTEPS = zeros(numberOfTemperatureDOFs,numInternalDOFs);
    KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
    %Initialize tangent for pre calculation of stresses
    Kinit = zeros(numInternalDOFs,numInternalDOFs);
    %Initialize sub-residual vector
    RT = zeros(numberOfTemperatureDOFs,1);
    RX = zeros(numberOfMechanicalDOFs,1);
    RSIG = zeros(numInternalDOFs,1);
    REPS = zeros(numInternalDOFs,1);
    HelmholtzN1 = 0;
    EntropyN1 = 0;
    InternalEnergyN1 = 0;    
    J = QREF(edofLocal,1:DIM)'*dNr';
    % Run through all Gauss points
    for k = 1:NGP
        indx = DIM*k-(DIM-1):DIM*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        dNx = (J(:,indx)')\dNr(indx,:);
        
        % Independent strains/conjungated stresses
        SIGMA_CN1v = zeros(6,1);
        CN1v = zeros(6,1);
        CNv = zeros(6,1);
        for m = 1:numNodes_C
            CNv = CNv + N_C(k,m)*CNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
            CN1v = CN1v + N_C(k,m)*CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            SIGMA_CN1v = SIGMA_CN1v + N_C(k,m)*SIGMA_CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        CN = CNv(map.VoigtInv);
        CN1 = CN1v(map.VoigtInv);
        SIGMA_CN1 = SIGMA_CN1v(map.VoigtInv);
        
        GN1v = zeros(6,1);
        SIGMA_GN1v=zeros(6,1);
        for m = 1:numNodes_G
            SIGMA_GN1v = SIGMA_GN1v + N_C(k,m)*SIGMA_GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            GN1v = GN1v + N_C(k,m)*GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        GN1 = GN1v(map.VoigtInv);
        SIGMA_GN1 = SIGMA_GN1v(map.VoigtInv);
        
        cN1 = 0;
        cN = 0;
        SIGMA_cN1 = 0;
        for m = 1:numNodes_c
            cN1 = cN1 + N_c(k,m)*cN1e(j,m);
            SIGMA_cN1 = SIGMA_cN1 + N_c(k,m)*SIGMA_cN1e(j,m);
            cN = cN + N_c(k,m)*cNe(j,m);
        end
        % Temperature
        thetaN1e = N(k,:)*thetaN1';
        thetaNe = N(k,:)*thetaN';
        
        % Deformation Gradient
        FxN1 = edN1*dNx';
        CxN1 = FxN1'*FxN1;
        
        % B-matrix (current configuration)
        BN1 = zeros(6,numberOfMechanicalDOFs);
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
        
        %% energy functions
        PsiIsoN1 = a*(trace(CN1) - 3) + b*(trace(GN1) - 3);
        PsiVolN1 = c/2*(sqrt(cN1) - 1)^2 - d*log(sqrt(cN1));
        PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
        PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*c*(sqrt(cN1) - 1);
        PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
        uThermN1 = kappa*(thetaN1e - theta0);
        uCoupledN1 = DIM*beta*theta0*c*(sqrt(cN1) - 1);
        etaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*c*(sqrt(cN1) - 1);
        uN1 = PsiIsoN1 + PsiVolN1 + uThermN1 + uCoupledN1;
        HelmholtzN1 = HelmholtzN1 + PsiN1*detJ*wp(k);
        EntropyN1 = EntropyN1 + etaN1*detJ*wp(k);
        InternalEnergyN1 = InternalEnergyN1 + uN1*detJ*wp(k);
        % derivatives
        Du_CN1 = a*eye(3);
        Du_GN1 = b*eye(3);
        Du_cN1 = c/2*(1 - cN1^(-1/2)) - d/(2*cN1) + DIM*beta*theta0*c/2*cN1^(-1/2);
        Du_thetaN1 = kappa;
        Deta_cN1 =  DIM*beta*c/2*cN1^(-1/2);
        Deta_thetaN1 = kappa/thetaN1e;
        D2u_cN1cN1 = c/4*cN1^(-3/2) + d/2*cN1^(-2) - DIM*beta*theta0*c/4*cN1^(-3/2);
        D2eta_cN1cN1 = -DIM*beta*c/4*cN1^(-3/2);
        D2eta_thetaN1thetaN1 = -kappa*thetaN1e^(-2);

        %% Residual
        % Balance of linear momentum
        RX = RX + BN1'*SIGMA_CN1v*detJ*wp(k);
        % Compatibility EPS
        rC = Du_CN1 - SIGMA_CN1 + wedge(SIGMA_GN1,CN1) + 1/3*SIGMA_cN1*GN1;
        rC = kron(N_C(k,:),Isym*rC(map.Voigt));
        rG = Du_GN1 - SIGMA_GN1 + 1/3*SIGMA_cN1*CN1;
        rG = kron(N_G(k,:),Isym*rG(map.Voigt));
        rc = kron(N_c(k,:),(Du_cN1 - thetaN1e*Deta_cN1 - SIGMA_cN1));
        REPS = REPS + [rC(:);rG(:);rc(:)]*detJ*wp(k);
        % Compatibility SIG
        rSIGC = CxN1 - CN1;
        rSIGC = kron(N_C(k,:),Isym*rSIGC(map.Voigt));
        rSIGG = 0.5*wedge(CN1,CN1)-GN1;
        rSIGG = kron(N_G(k,:),Isym*rSIGG(map.Voigt));
        rSIGc = kron(N_c(k,:),(1/3*sum(sum(GN1.*CN1))-cN1));
        RSIG = RSIG + [rSIGC(:);rSIGG(:);rSIGc(:)]*detJ*wp(k);
        % Energy Balance
        QHN1 = -k0*cN1^(-1)*GN1*(dNx*thetaN1');
        DotC = 1/DT*(CN1 - CN);
        rT = N(k,:)'*((thetaN1e-thetaNe)'/DT + (Deta_thetaN1)^(-1)*Deta_cN1*DotC(map.Voigt)'*(GN1(map.Voigt)'*Isym)') - Du_thetaN1^(-1)*dNx'*QHN1;
        RT = RT + rT*detJ*wp(k);

        %% Tangent K_X_o
        geometricalTangentOperator = 2*dNx'*SIGMA_CN1*dNx;
        kXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for g = 1:DIM
            kXX(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = geometricalTangentOperator;
        end
        KXX = KXX + kXX*detJ*wp(k);
        kXSIGC = kron(N_C(k,:),BN1');
        kXSIGG = zeros(numberOfMechanicalDOFs,6*numNodes_G);
        kXSIGc = zeros(numberOfMechanicalDOFs,numNodes_c);
        KXSIG = KXSIG + [kXSIGC kXSIGG kXSIGc]*detJ*wp(k);
        %% Tangent K_C_o
        kCSIGC = -kron(N_C(k,:)'*N_C(k,:),Isym);
        D = CN1;
        SecDiffOperator = [ 0          D(3,3)  D(2,2)      0        -D(3,2)  0;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1) 0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kCSIGG = kron(N_C(k,:)'*N_G(k,:),SecDiffOperator);
        D = SIGMA_GN1;
        SecDiffOperator = [ 0          D(3,3)  D(2,2)      0        -D(3,2)  0;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kCC = kron(N_C(k,:)'*N_C(k,:),SecDiffOperator);
        kCSIGc = 1/3*kron(N_C(k,:)'*N_c(k,:),Isym*GN1v);
        kCG = 1/3*SIGMA_cN1*kron(N_C(k,:)'*N_G(k,:),Isym);
        kCc = zeros(6*numNodes_C,numNodes_c);
        %% Tangent KG_o
        kGSIGG = -kron(N_G(k,:)'*N_G(k,:),Isym);
        kGSIGc = 1/3*kron(N_G(k,:)'*N_c(k,:),Isym*CN1v);
        kGC = 1/3*SIGMA_cN1*kron(N_G(k,:)'*N_C(k,:),Isym);
        kGSIGC = zeros(6*numNodes_G,6*numNodes_C);
        kGG = zeros(6*numNodes_G);
        kGc = zeros(6*numNodes_G,numNodes_c);
        %% Tangent K_c_o
        kcc = kron(N_c(k,:)'*N_c(k,:),1)*(D2u_cN1cN1 - thetaN1e*D2eta_cN1cN1);
        kcT = -N_c(k,:)'*Deta_cN1*N(k,:);
        kcSIGc = -kron(N_c(k,:)'*N_c(k,:),1);
        kcC = zeros(numNodes_c,6*numNodes_C);
        kcG = zeros(numNodes_c,6*numNodes_G);
        kcSIGC = zeros(numNodes_c,6*numNodes_C);
        kcSIGG = zeros(numNodes_c,6*numNodes_G);
        KEPSEPS = KEPSEPS + [kCC kCG kCc; kGC kGG kGc;kcC kcG kcc]*detJ*wp(k);
        KEPSSIG = KEPSSIG + [kCSIGC kCSIGG kCSIGc; kGSIGC kGSIGG kGSIGc; kcSIGC kcSIGG kcSIGc]*detJ*wp(k);
        KEPST = KEPST + [zeros(6*numNodes_C,numberOfTemperatureDOFs);zeros(6*numNodes_C,numberOfTemperatureDOFs);kcT]*detJ*wp(k);
        %% Tangent K_SIGC_o
        kSIGCX = kron(N_C(k,:)',BN1);
        kSIGCC = -kron(N_C(k,:)'*N_C(k,:),Isym);
        %% Tangent K_SIGG_o
        D = CN1;
        SecDiffOperator=[   0          D(3,3)  D(2,2)      0        -D(3,2)  0;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kSIGGC = kron(N_G(k,:)'*N_C(k,:),SecDiffOperator);
        kSIGGG = -kron(N_G(k,:)'*N_G(k,:),Isym);
        %% Tangent K_SIGc_o
        kSIGcC = 1/3*kron(N_C(k,:)'*N_c(k,:),(Isym*GN1v)');
        kSIGcG = 1/3*kron(N_G(k,:)'*N_c(k,:),(Isym*CN1v)');
        kSIGcc = -kron(N_c(k,:)'*N_c(k,:),1);
        kSIGGX = zeros(6*numNodes_G,numberOfMechanicalDOFs);
        kSIGcX = zeros(numNodes_c,numberOfMechanicalDOFs);
        kSIGCG = zeros(6*numNodes_C);
        kSIGCc = zeros(6*numNodes_C,numNodes_c);
        kSIGGc = zeros(6*numNodes_G,numNodes_c);
        KSIGX = KSIGX + [kSIGCX; kSIGGX; kSIGcX]*detJ*wp(k);
        KSIGEPS = KSIGEPS + [kSIGCC kSIGCG kSIGCc; kSIGGC kSIGGG kSIGGc;kSIGcC kSIGcG kSIGcc]*detJ*wp(k);
        %% Tangent K_T_o
        kTC = N(k,:)'*(Deta_thetaN1)^(-1)*Deta_cN1*kron(1/DT*N_C(k,:),(GN1(map.Voigt))'*Isym);
        mapdNx = [1;2;3;2;3;3];
        kTG = N(k,:)'*(Deta_thetaN1)^(-1)*Deta_cN1*kron(N_G(k,:),(DotC(map.Voigt))'*Isym) + Du_thetaN1^(-1)*k0*cN1^(-1)*kron((Isym*dNx(mapdNx,:))',N_G(k,:));
        kTc = N(k,:)'*(Deta_thetaN1)^(-1)*D2eta_cN1cN1*(DotC(map.Voigt)'*(GN1(map.Voigt)'*Isym)')*N_c(k,:) - k0*cN1^(-2)*(Du_thetaN1)^(-1)*(dNx'*(GN1*(dNx*thetaN1')))*N_c(k,:);
        DQHN1_thetaN1 = -k0*cN1^(-1)*GN1*dNx*eye(size(dNx,2));
        kTT = N(k,:)'*1/DT*N(k,:) - N(k,:)'*Deta_thetaN1^(-2)*(D2eta_thetaN1thetaN1*Deta_cN1)*N(k,:)*(sum(sum(DotC.*GN1))) - Du_thetaN1^(-1)*dNx'*DQHN1_thetaN1;
        KTEPS = KTEPS + [kTC kTG kTc]*detJ*wp(k);
        KTT = KTT + kTT*detJ*wp(k);

        ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1];
    end
    %% Static condensation
    ReCon1 = KEPSSIG\REPS;
    ReCon2 = KSIGEPS\RSIG;
    KeCon1 = KEPSSIG\KEPSEPS;
    KeCon2 = KSIGEPS\KSIGX;
    Re(mechanicalDOFs,1) = RX - KXSIG*ReCon1 + KXSIG*KeCon1*ReCon2;
    Re(temperatureDOFs,1) = RT - KTEPS*ReCon2;
    Ke(mechanicalDOFs,mechanicalDOFs) = KXX + KXSIG*KeCon1*KeCon2;
    Ke(mechanicalDOFs,temperatureDOFs) = -KXSIG*(KEPSSIG\KEPST);
    Ke(temperatureDOFs,mechanicalDOFs) = -KTEPS*(KSIGEPS\KSIGX);
    Ke(temperatureDOFs,temperatureDOFs) = KTT;
    %% Save values and matrices (for update)
    zw(j).RCon1 = ReCon1;
    zw(j).KCon1 = KeCon1;
    zw(j).RCon2 = ReCon2;
    zw(j).KCon2 = KeCon2;
    zw(j).KEPST = KEPST;
    zw(j).KEPSSIG = KEPSSIG;
    out(j).Re = Re;
    out(j).pK = Ke(:);
    out(j).ePot = ePot;
end
obj.StaConStruct = zw;
end
