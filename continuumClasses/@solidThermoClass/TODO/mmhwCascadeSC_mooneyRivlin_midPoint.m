function out = mmhwCascadeSC_mooneyRivlin_midPoint(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = mmhwCascadeSC_mooneyRivlinThermo_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% mmhwCascadeSC=: mixed method based with Hu-Washizu variational principle
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1/2 (midpoint scheme).
% 24.03.2017 M.S., A.J., M.F.

%% Check input
DT = 1;
TransRSig = 0;
LocalDofs = [];
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'TransRSig')
        TransRSig = 1;
    elseif strcmp(varargin{i},'LocalDofs')
        LocalDofs = varargin{i+1};
    end
end
ansatzfunction = obj.Ansatzfunction;
dNrAll = obj.SHAPEF.dNr;
%% Create residual and tangent
edof = obj.EDOF;
numElements = numel(edof);
if ~TransRSig
    edof = obj.GLOBALEDOF;
else
    edof = cell(1,numElements);
    for ii = 1:numElements
        zw = LocalDofs(obj.EDOF{ii},:)';
        edof{ii} = zw(:)';
    end
end
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
SIGMA_CNe  = obj.SIGMA_FN;
SIGMA_GNe  = obj.SIGMA_HN;
SIGMA_cNe  = obj.SIGMA_JN;
CNe  = obj.FN;
GNe  = obj.HN;
cNe  = obj.JNe;
% Element stuff
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];
Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
%% Shape functions
N = obj.SHAPEF.N;
wp = obj.SHAPEF.wp;
NGP = obj.NGP;
DIM = obj.DIM;
%% Material parameters
a = obj.MAT.c1;
b = obj.MAT.c2;
c = obj.MAT.c;
d = 2*(a+2*b);
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
theta0 = obj.MAT.theta0;
k0 = obj.MAT.k0;
DIM = obj.DIM;
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
%% Parfor
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
    dofTemp=4:DIM+1:numDOFs;
    dofMech=1:numDOFs;
    dofMech(4:DIM+1:numDOFs)=[];
    % Number of external dofs
    numTemp=numel(dofTemp);
    numMech=numel(dofMech);
     %Initialize element resdiual and tangent matrix
    Re=zeros(numDOFs,1);
    Ke=zeros(numDOFs);
    %Initialize sub-tangent matrices
    KXX = zeros(numMech,numMech);
    KXSIG = zeros(numMech,numInternalDOFs);
    KEPSEPS = zeros(numInternalDOFs,numInternalDOFs);
    KEPSSIG = zeros(numInternalDOFs,numInternalDOFs);
    KEPST = zeros(numInternalDOFs,numTemp);
    KSIGX = zeros(numInternalDOFs,numMech);
    KSIGEPS = zeros(numInternalDOFs,numInternalDOFs);
    KTEPS = zeros(numTemp,numInternalDOFs);
    KTT = zeros(numTemp,numTemp);
    %Initialize tangent for pre calculation of stresses
    Kinit = zeros(numInternalDOFs,numInternalDOFs);
    %Initialize sub-residual vector
    RT = zeros(numTemp,1);
    RX = zeros(numMech,1);
    RSIG = zeros(numInternalDOFs,1);
    REPS = zeros(numInternalDOFs,1);
    helmholtz = 0;
    entropyN1 = 0;
    internaleN1 = 0;
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
        SIGMA_CNv = zeros(6,1);
        CN1v = zeros(6,1);
        CNv = zeros(6,1);
        for m = 1:numNodes_C
            CN1v = CN1v + N_C(k,m)*CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            CNv = CNv + N_C(k,m)*CNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
            SIGMA_CN1v = SIGMA_CN1v + N_C(k,m)*SIGMA_CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            SIGMA_CNv = SIGMA_CNv + N_C(k,m)*SIGMA_CNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        CN1 = CN1v(map.VoigtInv);
        CN = CNv(map.VoigtInv);
        CN05 = 0.5*(CN1+CN);
        SIGMA_CN1 = SIGMA_CN1v(map.VoigtInv);
        SIGMA_CN = SIGMA_CNv(map.VoigtInv);
        SIGMA_CN05 = 0.5*(SIGMA_CN1+SIGMA_CN);
        GN1v = zeros(6,1);
        GNv = zeros(6,1);
        SIGMA_GN1v=zeros(6,1);
        SIGMA_GNv=zeros(6,1);
        for m = 1:numNodes_G
            SIGMA_GN1v = SIGMA_GN1v + N_C(k,m)*SIGMA_GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            SIGMA_GNv = SIGMA_GNv + N_C(k,m)*SIGMA_GNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
            GN1v = GN1v + N_C(k,m)*GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            GNv = GNv + N_C(k,m)*GNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        GN1 = GN1v(map.VoigtInv);
        GN = GNv(map.VoigtInv);
        GN05 = 0.5*(GN1+GN);
        SIGMA_GN1 = SIGMA_GN1v(map.VoigtInv);
        SIGMA_GN = SIGMA_GNv(map.VoigtInv);
        SIGMA_GN05 = 0.5*(SIGMA_GN1+SIGMA_GN);
        cN1 = 0;
        cN = 0;
        SIGMA_cN1 = 0;
        SIGMA_cN = 0;
        for m = 1:numNodes_c
            cN1 = cN1 + N_c(k,m)*cN1e(j,m);
            SIGMA_cN1 = SIGMA_cN1 + N_c(k,m)*SIGMA_cN1e(j,m);
            SIGMA_cN = SIGMA_cN + N_c(k,m)*SIGMA_cNe(j,m);
            cN = cN + N_c(k,m)*cNe(j,m);
        end
        cN05 = 0.5*(cN1+cN);
        SIGMA_cN05 = 0.5*(SIGMA_cN1+SIGMA_cN);
        % Temperature
        thetaN1e = N(k,:)*thetaN1';
        thetaNe = N(k,:)*thetaN';
        thetaN05  = 0.5*(thetaN1+thetaN);
        thetaN05e = N(k,:)*thetaN05';
        edN05 = 0.5*(edN1+edN);
        % Deformation Gradient
        FxN05 = edN05*dNx';
        CxN05 = FxN05'*FxN05;
        % B-matrix (current configuration)
        BN05 = zeros(6,numMech);
        BN05(1,1:3:end) = FxN05(1,1)*dNx(1,:);
        BN05(1,2:3:end) = FxN05(2,1)*dNx(1,:);
        BN05(1,3:3:end) = FxN05(3,1)*dNx(1,:);
        BN05(2,1:3:end) = FxN05(1,2)*dNx(2,:);
        BN05(2,2:3:end) = FxN05(2,2)*dNx(2,:);
        BN05(2,3:3:end) = FxN05(3,2)*dNx(2,:);
        BN05(3,1:3:end) = FxN05(1,3)*dNx(3,:);
        BN05(3,2:3:end) = FxN05(2,3)*dNx(3,:);
        BN05(3,3:3:end) = FxN05(3,3)*dNx(3,:);
        BN05(4,1:3:end) = FxN05(1,1)*dNx(2,:)+FxN05(1,2)*dNx(1,:);
        BN05(4,2:3:end) = FxN05(2,1)*dNx(2,:)+FxN05(2,2)*dNx(1,:);
        BN05(4,3:3:end) = FxN05(3,1)*dNx(2,:)+FxN05(3,2)*dNx(1,:);
        BN05(5,1:3:end) = FxN05(1,2)*dNx(3,:)+FxN05(1,3)*dNx(2,:);
        BN05(5,2:3:end) = FxN05(2,2)*dNx(3,:)+FxN05(2,3)*dNx(2,:);
        BN05(5,3:3:end) = FxN05(3,2)*dNx(3,:)+FxN05(3,3)*dNx(2,:);
        BN05(6,1:3:end) = FxN05(1,1)*dNx(3,:)+FxN05(1,3)*dNx(1,:);
        BN05(6,2:3:end) = FxN05(2,1)*dNx(3,:)+FxN05(2,3)*dNx(1,:);
        BN05(6,3:3:end) = FxN05(3,1)*dNx(3,:)+FxN05(3,3)*dNx(1,:);
        BN05 = 2*BN05;
        % Strain energy functions
        PsiIsoN1 = a*(trace(CN1)-3) + b*(trace(GN1)-3);
        PsiVolN1 = -d*log(sqrt(cN1)) + c/2*(sqrt(cN1)-1)^2;
        PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
        PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c*(sqrt(cN1) - 1) - d/sqrt(cN1));
        PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
        emechN1 = PsiVolN1;
        emixN1 = kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c*(sqrt(cN1) - 1) - d/sqrt(cN1));
        etotalN1 = PsiIsoN1 + emechN1 + emixN1;
        etaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*(c*(sqrt(cN1)-1)-d/sqrt(cN1));
        helmholtz = helmholtz + PsiN1*detJ*wp(k);
        entropyN1 = entropyN1 + etaN1*detJ*wp(k);
        internaleN1 = internaleN1 + etotalN1*detJ*wp(k);
        %% Impuls balance
        % Derivative of the strain energy function
        DPsi_C = a*eye(3);
        DPsi_G = b*eye(3);
        DPsi_cN05 = -d/(2*cN05) + c/2*(1-1/sqrt(cN05)) - DIM*beta*(thetaN05e - theta0)*(c/2*1/sqrt(cN05) + d/2*(cN05)^(-3/2));
        % Hessian operator of the complementary energy function
        D2Psi_cc = d/(2*cN05^2) + c/(4*(cN05)^(3/2)) - DIM*beta*(thetaN05e - theta0)*(-c/4*cN05^(-3/2) - 3*d/4*cN05^(-5/2));
        D2Psi_ct = 0.5*DIM*beta*(c/2*1/sqrt(cN05)+d/2*(cN05)^(-3/2));
        %rX
        rX = BN05'*SIGMA_CN05(map.Voigt);
        %% Compatibility EPS
        %rEPS
        DPsi_C_rC = [DPsi_C(1,1) DPsi_C(2,2) DPsi_C(3,3) 2*DPsi_C(1,2) 2*DPsi_C(2,3) 2*DPsi_C(1,3)].';
        SIGMA_CN05_rC = [SIGMA_CN05(1,1) SIGMA_CN05(2,2) SIGMA_CN05(3,3) 2*SIGMA_CN05(1,2) 2*SIGMA_CN05(2,3) 2*SIGMA_CN05(1,3)].';
        wedge_SIGMA_GN05CN05 = wedge(SIGMA_GN05,CN05);
        wedge_SIGMA_GN05CN05_rC = [wedge_SIGMA_GN05CN05(1,1) wedge_SIGMA_GN05CN05(2,2) wedge_SIGMA_GN05CN05(3,3) 2*wedge_SIGMA_GN05CN05(1,2) 2*wedge_SIGMA_GN05CN05(2,3) 2*wedge_SIGMA_GN05CN05(1,3)].';
        GN05_rC = [GN05(1,1) GN05(2,2) GN05(3,3) 2*GN05(1,2) 2*GN05(2,3) 2*GN05(1,3)].';
        DPsi_G_rG = [DPsi_G(1,1) DPsi_G(2,2) DPsi_G(3,3) 2*DPsi_G(1,2) 2*DPsi_G(2,3) 2*DPsi_G(1,3)].';
        SIGMA_GN05_rG = [SIGMA_GN05(1,1) SIGMA_GN05(2,2) SIGMA_GN05(3,3) 2*SIGMA_GN05(1,2) 2*SIGMA_GN05(2,3) 2*SIGMA_GN05(1,3)].';
        CN05_rG = [CN05(1,1) CN05(2,2) CN05(3,3) 2*CN05(1,2) 2*CN05(2,3) 2*CN05(1,3)].';
        if TransRSig==1
            rI = kron(N_c(k,:),(DPsi_cN05));
            rI = rI(:);
            rC = kron(N_C(k,:),(DPsi_C_rC+wedge_SIGMA_GN05CN05_rC+1/3*SIGMA_cN05*GN05_rC));
            rC = rC(:);
            rG = kron(N_G(k,:),(DPsi_G_rG+1/3*SIGMA_cN05*CN05_rG));
            rG = rG(:);
        else
            rI = kron(N_c(k,:),(DPsi_cN05-SIGMA_cN05));
            rI = rI(:);
            rC = kron(N_C(k,:),(DPsi_C_rC-SIGMA_CN05_rC+wedge_SIGMA_GN05CN05_rC+1/3*SIGMA_cN05*GN05_rC));
            rC = rC(:);
            rG = kron(N_G(k,:),(DPsi_G_rG-SIGMA_GN05_rG+1/3*SIGMA_cN05*CN05_rG));
            rG = rG(:);
        end
        %% Compatibility SIG
        %rSIG
        CxN05_rSC = [CxN05(1,1) CxN05(2,2) CxN05(3,3) 2*CxN05(1,2) 2*CxN05(2,3) 2*CxN05(1,3)].';
        CN05_rSC = [CN05(1,1) CN05(2,2) CN05(3,3) 2*CN05(1,2) 2*CN05(2,3) 2*CN05(1,3)].';
        wedge_CN05CN05 = wedge(CN05,CN05);
        wedge_CN05CN05_rSG = [wedge_CN05CN05(1,1) wedge_CN05CN05(2,2) wedge_CN05CN05(3,3) 2*wedge_CN05CN05(1,2) 2*wedge_CN05CN05(2,3) 2*wedge_CN05CN05(1,3)].';
        GN05_rSG =  [GN05(1,1) GN05(2,2) GN05(3,3) 2*GN05(1,2) 2*GN05(2,3) 2*GN05(1,3)].';
        
        rSIGC = kron(N_C(k,:),CxN05_rSC-CN05_rSC);
        rSIGC = rSIGC(:);
        rSIGG = kron(N_G(k,:),(0.5*wedge_CN05CN05_rSG-GN05_rSG));
        rSIGG = rSIGG(:);
        rSIGI = kron(N_c(k,:),(1/3*sum(sum(GN05.*CN05))-cN05));
        rSIGI = rSIGI(:);
        %% Energy Balance
        % Entropy
        etaN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c*(sqrt(cN1)-1)-d/sqrt(cN1));
        etaN = kappa*log(thetaNe/theta0)+DIM*beta*(c*(sqrt(cN)-1)-d/sqrt(cN));
        % Fourier law
        Q = -k0*(1/cN05*GN05*(dNx*thetaN05'));
        % rT
        rT = (N(k,:)'*thetaN05e*(etaN1-etaN)/DT-dNx'*Q);
        %% Tangent K_X_o
        kXX = zeros(numMech,numMech);
        A1 = 0.5*2*dNx'*SIGMA_CN05*dNx;
        for g = 1:DIM
            kXX(g:DIM:numMech,g:DIM:numMech) = A1;
        end
        kXSIGC = 0.5*kron(N_C(k,:),BN05');
        %% Tangent K_C_o
        kCSIGC = -0.5*kron(N_C(k,:)'*N_C(k,:),Isym);
        D = CN05;
        SecDiffOperator = [ 0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1) 0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kCSIGG = 0.5*kron(N_C(k,:)'*N_G(k,:),SecDiffOperator);
        D = SIGMA_GN05;
        SecDiffOperator = [ 0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kCC = 0.5*kron(N_C(k,:)'*N_C(k,:),SecDiffOperator);
        kCSIGc = 0.5*1/3*kron(N_C(k,:)'*N_c(k,:),GN05_rC);
        kCG = 0.5*1/3*SIGMA_cN05*kron(N_C(k,:)'*N_G(k,:),Isym);
        %% Tangent KG_o
        kGSIGG = -0.5*kron(N_G(k,:)'*N_G(k,:),Isym);
        kGSIGc = 0.5*1/3*kron(N_G(k,:)'*N_c(k,:),CN05_rG);
        kGC = 0.5*1/3*SIGMA_cN05*kron(N_G(k,:)'*N_C(k,:),Isym);
        %% Tangent K_c_o
        kcc = 0.5*kron(N_c(k,:)'*N_c(k,:),1)*D2Psi_cc;
        kcSIGc = -0.5*kron(N_c(k,:)'*N_c(k,:),1);
        %% Tangent K_SIGC_o
        kSIGCX = 0.5*kron(N_C(k,:)',BN05);
        kSIGCC = -0.5*kron(N_C(k,:)'*N_C(k,:),Isym);
        %% Tangent K_SIGG_o
        D=CN05;
        SecDiffOperator=[   0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kSIGGC = 0.5*kron(N_G(k,:)'*N_C(k,:),SecDiffOperator);
        kSIGGG = -0.5*kron(N_G(k,:)'*N_G(k,:),Isym);
        %% Tangent K_SIGc_o
        kSIGcG = 0.5*1/3*kron(N_G(k,:),N_c(k,:)'*CN05_rG');
        kSIGcC = 0.5*1/3*kron(N_C(k,:),N_c(k,:)'*GN05_rC');
        kSIGcc = -0.5*kron(N_c(k,:)'*N_c(k,:),1);
        %% Tangent K_TT
        eta_thetaN1 = kappa/thetaN1e;
        Q_The = -0.5*k0*(1/cN05*GN05*dNx*eye(size(dNx,2)));
        kTT =  N(k,:)'*0.5*(etaN1-etaN)*N(k,:)/DT+thetaN05e*N(k,:)'*eta_thetaN1*N(k,:)/DT-dNx'*Q_The;
        %% Tangent K_TG
        I = eye(3);
        Q_1G=zeros(3,3);
        Q_2G=zeros(3,3);
        Q_3G=zeros(3,3);
        Q_1G(:,:)=-0.5*k0*1/cN05*(I(:,1)*(dNx*thetaN05')'+(dNx*thetaN05')*I(1,:));
        Q_2G(:,:)=-0.5*k0*1/cN05*(I(:,2)*(dNx*thetaN05')'+(dNx*thetaN05')*I(2,:));
        Q_3G(:,:)=-0.5*k0*1/cN05*(I(:,3)*(dNx*thetaN05')'+(dNx*thetaN05')*I(3,:));
        Q_Gv = zeros(DIM,6);
        Q_Gv(1,:)=Q_1G(map.Voigt);
        Q_Gv(2,:)=Q_2G(map.Voigt);
        Q_Gv(3,:)=Q_3G(map.Voigt);
        kTG = kron(N_G(k,:),-dNx'*0.5*Q_Gv*Isym);
        %% Tangent K_Tc
        Q_c = 0.5*k0*(1/cN05^2*GN05*(dNx*thetaN05'));
        eta_cN1 = DIM*beta*(c/2*1/sqrt(cN1)+d/2*cN1^(-3/2));
        kTc = N(k,:)'*thetaN05e*eta_cN1*N_c(k,:)/DT-dNx'*Q_c*N_c(k,:);
        kcT = -N_c(k,:)'*D2Psi_ct*N(k,:);
        %% Zeros
        kXSIGG = zeros(numMech,6*numNodes_G);
        kXSIGc = zeros(numMech,numNodes_c);
        kCc = zeros(6*numNodes_C,numNodes_c);
        kGG = zeros(6*numNodes_G);
        kGc = zeros(6*numNodes_G,numNodes_c);
        kcC = zeros(numNodes_c,6*numNodes_C);
        kcG = zeros(numNodes_c,6*numNodes_G);
        kGSIGC = zeros(6*numNodes_G,6*numNodes_C);
        kcSIGC = zeros(numNodes_c,6*numNodes_C);
        kcSIGG = zeros(numNodes_c,6*numNodes_G);
        kSIGGX = zeros(6*numNodes_G,numMech);
        kSIGcX = zeros(numNodes_c,numMech);
        kSIGCG = zeros(6*numNodes_C);
        kSIGCc = zeros(6*numNodes_C,numNodes_c);
        kSIGGc = zeros(6*numNodes_G,numNodes_c);
        % Tangents for pre calculation of stresses
        kSIGCSIGC = kron(N_C(k,:)'*N_C(k,:),Isym);
        kSIGGSIGG = kron(N_G(k,:)'*N_G(k,:),Isym);
        kSIGcSIGc = N_c(k,:)'*N_c(k,:);
        %% summation of residual
        RX = RX + rX*detJ*wp(k);
        REPS = REPS + [rC;rG;rI]*detJ*wp(k);
        RSIG = RSIG + [rSIGC;rSIGG;rSIGI]*detJ*wp(k);
        RT = RT + rT*detJ*wp(k);
        %% summation of tangents
        KXX = KXX + kXX*detJ*wp(k);
        KXSIG = KXSIG + [kXSIGC kXSIGG kXSIGc]*detJ*wp(k);
        KEPSEPS = KEPSEPS + [kCC kCG kCc; kGC kGG kGc;kcC kcG kcc]*detJ*wp(k);
        KEPSSIG = KEPSSIG + [kCSIGC kCSIGG kCSIGc; kGSIGC kGSIGG kGSIGc; kcSIGC kcSIGG kcSIGc]*detJ*wp(k);
        KSIGX = KSIGX + [kSIGCX; kSIGGX; kSIGcX]*detJ*wp(k);
        KSIGEPS = KSIGEPS + [kSIGCC kSIGCG kSIGCc; kSIGGC kSIGGG kSIGGc;kSIGcC kSIGcG kSIGcc]*detJ*wp(k);
        KTT = KTT + kTT*detJ*wp(k);
        KTEPS = KTEPS + [zeros(numTemp,6*numNodes_C) kTG kTc]*detJ*wp(k);
        KEPST = KEPST + [zeros(6*numNodes_C,numTemp);zeros(6*numNodes_C,numTemp);kcT]*detJ*wp(k);
        % Tangents for pre calculation of stresses
        if TransRSig == 1
            KSIGCSIGC = kSIGCSIGC*detJ*wp(k);
            KSIGGSIGG = kSIGGSIGG*detJ*wp(k);
            KSIGcSIGc = kSIGcSIGc*detJ*wp(k);
            Kinit = Kinit + [KSIGCSIGC zeros(6*numNodes_C,6*numNodes_C) zeros(6*numNodes_C,numNodes_c);zeros(6*numNodes_G,6*numNodes_C) KSIGGSIGG zeros(6*numNodes_G,numNodes_c);zeros(numNodes_c,6*numNodes_C) zeros(numNodes_c,6*numNodes_G) KSIGcSIGc];
        end
    end
    if TransRSig==1
        SIG = Kinit\REPS;
        out(j).Re = SIG;
        out(j).pK = Ke(:);
        out(j).ePot = 0;
    else
        %% Static condesation process
        ReCon1 = KEPSSIG\REPS;
        ReCon2 = KSIGEPS\RSIG;
        KeCon1 = KEPSSIG\KEPSEPS;
        KeCon2 = KSIGEPS\KSIGX;
        Re(dofMech,1)=RX-KXSIG*ReCon1 + KXSIG*KeCon1*ReCon2;
        Re(dofTemp,1)=RT-KTEPS*ReCon2;
        Ke(dofMech,dofMech)=KXX+KXSIG*KeCon1*KeCon2;
        Ke(dofMech,dofTemp)=-KXSIG*(KEPSSIG\KEPST);
        Ke(dofTemp,dofMech)=-KTEPS*(KSIGEPS\KSIGX);
        Ke(dofTemp,dofTemp)=KTT;
        %% Save values and matrices (for update)
        zw(j).RCon1 = ReCon1;
        zw(j).KCon1 = KeCon1;
        zw(j).RCon2 = ReCon2;
        zw(j).KCon2 = KeCon2;
        zw(j).KEPST = KEPST;
        zw(j).KEPSSIG = KEPSSIG;
        out(j).Re = Re;
        out(j).pK = Ke(:);
        out(j).ePot = [helmholtz entropyN1 internaleN1];
    end
end
obj.StaConStruct = zw;
end
