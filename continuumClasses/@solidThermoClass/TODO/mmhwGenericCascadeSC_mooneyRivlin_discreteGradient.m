function out = mmhwGenericCascadeSC_mooneyRivlin_discreteGradient(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = mmhwGenericCascadeSC_mooneyRivlinThermo_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% mmhwCascadeSC=: mixed method based with Hu-Washizu variational principle
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% discreteGradient.
% 09.06.2017 M.S., A.J., M.F.

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
GNe  = obj.HN;
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
    edN05 = 0.5*(edN + edN1);
    thetaN1 = QN1(edofLocal,DIM+1)';
    thetaN = QN(edofLocal,DIM+1)';
    thetaN05 = 0.5*(thetaN + thetaN1);
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
    DeltaU = 0;
    ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 DeltaU];  
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
        CN05v=zeros(6,1);
        for m = 1:numNodes_C
            CNv = CNv + N_C(k,m)*CNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
            CN1v = CN1v + N_C(k,m)*CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            CN05v= CN05v + N_C(k,m)*0.5*(CN1e(j,ones(1,6)*(m-1)*6+(1:6))'+CNe(j,ones(1,6)*(m-1)*6+(1:6))');
            SIGMA_CN1v = SIGMA_CN1v + N_C(k,m)*SIGMA_CN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        CN = CNv(map.VoigtInv);
        CN1 = CN1v(map.VoigtInv);
        CN05 = CN05v(map.VoigtInv);
        SIGMA_CN1 = SIGMA_CN1v(map.VoigtInv);
        
        GNv = zeros(6,1);
        GN1v = zeros(6,1);
        GN05v=zeros(6,1);
        SIGMA_GN1v=zeros(6,1);
        for m = 1:numNodes_G
            GNv = GNv + N_G(k,m)*GNe(j,ones(1,6)*(m - 1)*6 + (1:6))';
            GN1v = GN1v + N_G(k,m)*GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
            GN05v= GN05v + N_G(k,m)*0.5*(GN1e(j,ones(1,6)*(m-1)*6+(1:6))'+GNe(j,ones(1,6)*(m-1)*6+(1:6))');
            SIGMA_GN1v = SIGMA_GN1v + N_G(k,m)*SIGMA_GN1e(j,ones(1,6)*(m - 1)*6 + (1:6))';
        end
        GN = GNv(map.VoigtInv);        
        GN1 = GN1v(map.VoigtInv);
        GN05 = GN05v(map.VoigtInv);
        SIGMA_GN1 = SIGMA_GN1v(map.VoigtInv);
        
        cN1 = 0;
        cN = 0;
        cN05 = 0;
        SIGMA_cN1 = 0;
        for m = 1:numNodes_c
            cN1 = cN1 + N_c(k,m)*cN1e(j,m);
            cN = cN + N_c(k,m)*cNe(j,m);
            cN05 = cN05 + N_c(k,m)*0.5*(cN1e(j,m)+cNe(j,m));            
            SIGMA_cN1 = SIGMA_cN1 + N_c(k,m)*SIGMA_cN1e(j,m);
        end
        % Temperature
        thetaN1e = N(k,:)*thetaN1';
        thetaNe = N(k,:)*thetaN';
        thetaN05e = N(k,:)*thetaN05';
        
        % Deformation Gradient
        FxN1 = edN1*dNx';
        FxN05 = edN05*dNx';
        CxN1 = FxN1'*FxN1;
        Galgo =1/3*(wedge(CN05,CN05) + GN05);

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
    
        % B-matrix (current configuration)
        BN05 = zeros(6,numberOfMechanicalDOFs);
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
        
        %% energy functions
        PsiIsoN1 = a*(trace(CN1) - 3) + b*(trace(GN1) - 3);
        PsiIsoN = a*(trace(CN) - 3) + b*(trace(GN) - 3);
        PsiVolN1 = c/2*(sqrt(cN1) - 1)^2 - d*log(sqrt(cN1));
        PsiVolN = c/2*(sqrt(cN) - 1)^2 - d*log(sqrt(cN));
        PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
        PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*c*(sqrt(cN1) - 1);
        PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
        uThermN1 = kappa*(thetaN1e - theta0);
        uThermN = kappa*(thetaNe - theta0);
        uCoupledN1 = DIM*beta*theta0*c*(sqrt(cN1) - 1);
        uCoupledN = DIM*beta*theta0*c*(sqrt(cN) - 1);
        etaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*c*(sqrt(cN1) - 1);        
        uN1 = PsiIsoN1 + PsiVolN1 + uThermN1 + uCoupledN1;
        uN = PsiIsoN + PsiVolN + uThermN + uCoupledN;
        HelmholtzN1 = HelmholtzN1 + PsiN1*detJ*wp(k);
        EntropyN1 = EntropyN1 + etaN1*detJ*wp(k);
        InternalEnergyN1 = InternalEnergyN1 + uN1*detJ*wp(k);
        DeltaU = DeltaU + (uN1 - uN)*detJ*wp(k);
        
        % derivatives
%         Du_CN1 = a*eye(3);
%         Du_GN1 = b*eye(3);
%         Du_cN1 = c/2*(1 - cN1^(-1/2)) - d/(2*cN1) + DIM*beta*theta0*c/2*cN1^(-1/2);
%         Du_thetaN1 = kappa;
%         Deta_cN1 =  DIM*beta*c/2*cN1^(-1/2);
%         Deta_thetaN1 = kappa/thetaN1e;
%         D2u_cN1cN1 = c/4*cN1^(-3/2) + d/2*cN1^(-2) - DIM*beta*theta0*c/4*cN1^(-3/2);
%         D2eta_cN1cN1 = -DIM*beta*c/4*cN1^(-3/2);
%         D2eta_thetaN1thetaN1 = -kappa*thetaN1e^(-2);
        Du_CN05 = a*eye(3);
        Du_GN05 = b*eye(3);
        Du_theta = kappa;

        %% Residual
        % Balance of linear momentum
        RX = RX + BN05'*SIGMA_CN1v*detJ*wp(k);
        % Compatibility EPS
        etacN1thetaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*c*(sqrt(cN1)-1);
        etacN1thetaN = kappa*log(thetaNe/theta0) + DIM*beta*c*(sqrt(cN1)-1);
        etacNthetaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*c*(sqrt(cN)-1);
        etacNthetaN = kappa*log(thetaNe/theta0) + DIM*beta*(c*(sqrt(cN)-1));
        ucN1thetaN1 =   a*(trace(CN1) - 3) + b*(trace(GN1) - 3) + c/2*(sqrt(cN1) - 1)^2 - d*log(sqrt(cN1)) + kappa*(thetaN1e - theta0) + DIM*beta*theta0*c*(sqrt(cN1) - 1);
        ucN1thetaN =    a*(trace(CN1) - 3) + b*(trace(GN1) - 3) + c/2*(sqrt(cN1) - 1)^2 - d*log(sqrt(cN1)) + kappa*(thetaNe - theta0) + DIM*beta*theta0*c*(sqrt(cN1) - 1);
        ucNthetaN1 =    a*(trace(CN1) - 3) + b*(trace(GN1) - 3) + c/2*(sqrt(cN) - 1)^2 - d*log(sqrt(cN)) + kappa*(thetaN1e - theta0) + DIM*beta*theta0*c*(sqrt(cN) - 1);
        ucNthetaN =     a*(trace(CN1) - 3) + b*(trace(GN1) - 3) + c/2*(sqrt(cN) - 1)^2 - d*log(sqrt(cN)) + kappa*(thetaNe - theta0) + DIM*beta*theta0*c*(sqrt(cN) - 1);        
        DucN1thetaN_cN1 = (c*(cN1^(1/2) - 1))/(2*cN1^(1/2)) - d/(2*cN1) + DIM*beta*theta0*(c/(2*cN1^(1/2)));
        DucN1thetaN1_cN1 = (c*(cN1^(1/2) - 1))/(2*cN1^(1/2)) - d/(2*cN1) + DIM*beta*theta0*(c/(2*cN1^(1/2)));
        DetacN1thetaN_cN1 =  DIM*beta*(c/(2*cN1^(1/2)));
        DetacN1thetaN1_cN1 =  DIM*beta*(c/(2*cN1^(1/2)));
        DetacN1thetaN1_theta = kappa/thetaN1e;
        DetacNthetaN1_theta = kappa/thetaN1e;
        if abs(thetaN1e - thetaNe) >= 10*eps
            Deta_thetaLeft = (etacN1thetaN1 - etacN1thetaN )*(thetaN1e - thetaNe)^(-1);
            Deta_thetaRight = (etacNthetaN1 - etacNthetaN )*(thetaN1e - thetaNe)^(-1);
            Deta_theta = 1/2*(Deta_thetaLeft + Deta_thetaRight);
            D2eta_thetatheta = 1/2*((DetacN1thetaN1_theta + DetacNthetaN1_theta)/(thetaN1e - thetaNe) + (-(thetaN1e - thetaNe)^(-2))*(etacN1thetaN1 - etacN1thetaN + etacNthetaN1 - etacNthetaN));
            D2eta_thetatheta = 2*D2eta_thetatheta;
        else
            Deta_theta = kappa/thetaN05e;
            D2eta_thetatheta = -kappa*thetaN05e^(-2);
        end
        thetaAlgo = Deta_theta^(-1)*Du_theta;
        DthetaAlgo_theta = -Deta_theta^(-2)*D2eta_thetatheta*Du_theta;
        if abs(cN1 - cN) >= 10*eps
            Du_cLeft = (ucN1thetaN-ucNthetaN)/(cN1 - cN);
            Du_cRight = (ucN1thetaN1-ucNthetaN1)/(cN1 - cN);
            Du_c = 0.5*(Du_cLeft + Du_cRight);
            D2u_cc = ((DucN1thetaN_cN1 + DucN1thetaN1_cN1)/(cN1 - cN) + (-(cN1 - cN)^(-2))*(ucN1thetaN - ucNthetaN + ucN1thetaN1 - ucNthetaN1));
            Deta_cLeft = (etacN1thetaN-etacNthetaN)/(cN1 - cN);
            Deta_cRight = (etacN1thetaN1-etacNthetaN1)/(cN1 - cN);
            Deta_c = 0.5*(Deta_cLeft+Deta_cRight);
            D2eta_cc = ((DetacN1thetaN_cN1 + DetacN1thetaN1_cN1)/(cN1 - cN) + (-(cN1 - cN)^(-2))*(etacN1thetaN - etacNthetaN + etacN1thetaN1 - etacNthetaN1));
        else
            Du_c = c/2*(1 - 1/sqrt(cN05)) - d/(2*cN05) + DIM*beta*theta0*(c/2*1/sqrt(cN05));
            D2u_cc = c/(4*cN05^(3/2)) + d/(2*cN05^2) - DIM*beta*theta0*(c/(4*cN05^(3/2)));
            Deta_c =  DIM*beta*(c/2*1/sqrt(cN05));        
            D2eta_cc = -DIM*beta*(c/(4*cN05^(3/2)));
        end
        
        rC = Du_CN05 - SIGMA_CN1 + wedge(SIGMA_GN1,CN05) + 1/3*SIGMA_cN1*GN05;
        rC = kron(N_C(k,:),Isym*rC(map.Voigt));
        rG = Du_GN05 - SIGMA_GN1 + 1/3*SIGMA_cN1*CN05;
        rG = kron(N_G(k,:),Isym*rG(map.Voigt));
        rc = kron(N_c(k,:),(Du_c - thetaAlgo*Deta_c - SIGMA_cN1));
        REPS = REPS + [rC(:);rG(:);rc(:)]*detJ*wp(k);
        % Compatibility SIG
        rSIGC = CxN1 - CN1;
        rSIGC = kron(N_C(k,:),Isym*rSIGC(map.Voigt));
        rSIGG = 0.5*wedge(CN1,CN1)-GN1;
        rSIGG = kron(N_G(k,:),Isym*rSIGG(map.Voigt));
        rSIGc = kron(N_c(k,:),(1/3*sum(sum(GN1.*CN1))-cN1));
        RSIG = RSIG + [rSIGC(:);rSIGG(:);rSIGc(:)]*detJ*wp(k);
        % Energy Balance
        QHN05 = -k0*cN05^(-1)*GN05*(dNx*thetaN05');
        DotC = 1/DT*(CN1 - CN);
        rT = N(k,:)'*((thetaN1e-thetaNe)'/DT + (Deta_theta)^(-1)*Deta_c*DotC(map.Voigt)'*(Galgo(map.Voigt)'*Isym)') - Du_theta^(-1)*dNx'*QHN05;
        RT = RT + rT*detJ*wp(k);

        %% Tangent K_X_o
        geometricalTangentOperator = 2*dNx'*SIGMA_CN1*dNx;
        kXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for g = 1:DIM
            kXX(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = geometricalTangentOperator;
        end
        KXX = KXX + 0.5*kXX*detJ*wp(k);
        kXSIGC = kron(N_C(k,:),BN05');
        kXSIGG = zeros(numberOfMechanicalDOFs,6*numNodes_G);
        kXSIGc = zeros(numberOfMechanicalDOFs,numNodes_c);
        KXSIG = KXSIG + [kXSIGC kXSIGG kXSIGc]*detJ*wp(k);
        %% Tangent K_C_o
        kCSIGC = -kron(N_C(k,:)'*N_C(k,:),Isym);
        D = CN05;
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
        kCSIGc = 1/3*kron(N_C(k,:)'*N_c(k,:),Isym*GN05v);
        kCG = 1/3*SIGMA_cN1*kron(N_C(k,:)'*N_G(k,:),Isym);
        kCc = zeros(6*numNodes_C,numNodes_c);
        %% Tangent KG_o
        kGSIGG = -kron(N_G(k,:)'*N_G(k,:),Isym);
        kGSIGc = 1/3*kron(N_G(k,:)'*N_c(k,:),Isym*CN05v);
        kGC = 1/3*SIGMA_cN1*kron(N_G(k,:)'*N_C(k,:),Isym);
        kGSIGC = zeros(6*numNodes_G,6*numNodes_C);
        kGG = zeros(6*numNodes_G);
        kGc = zeros(6*numNodes_G,numNodes_c);
        %% Tangent K_c_o
        kcc = (kron(N_c(k,:)'*N_c(k,:),1)*(D2u_cc - thetaAlgo*D2eta_cc));
        kcT = -N_c(k,:)'*DthetaAlgo_theta*Deta_c*N(k,:);
        kcSIGc = -kron(N_c(k,:)'*N_c(k,:),1);
        kcC = zeros(numNodes_c,6*numNodes_C);
        kcG = zeros(numNodes_c,6*numNodes_G);
        kcSIGC = zeros(numNodes_c,6*numNodes_C);
        kcSIGG = zeros(numNodes_c,6*numNodes_G);
        KEPSEPS = KEPSEPS + 0.5*[kCC kCG kCc; kGC kGG kGc;kcC kcG kcc]*detJ*wp(k);
        KEPSSIG = KEPSSIG + [kCSIGC kCSIGG kCSIGc; kGSIGC kGSIGG kGSIGc; kcSIGC kcSIGG kcSIGc]*detJ*wp(k);
        KEPST = KEPST + 0.5*[zeros(6*numNodes_C,numberOfTemperatureDOFs);zeros(6*numNodes_C,numberOfTemperatureDOFs);kcT]*detJ*wp(k);
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
        D = CN05;
        SecDiffOperator=[   0          D(3,3)  D(2,2)      0        -D(3,2)  0;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)    0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        kTC = N(k,:)'*(Deta_theta)^(-1)*Deta_c*kron(1/DT*N_C(k,:),(Galgo(map.Voigt))'*Isym) + (Deta_theta)^(-1)*Deta_c*1/3*kron(N(k,:)'*N_C(k,:),DotC(map.Voigt)'*SecDiffOperator);
        mapdNx = [1;2;3;2;3;3];
        kTG = N(k,:)'*(Deta_theta)^(-1)*Deta_c*kron(N_G(k,:),(DotC(map.Voigt))'*1/3*Isym) + Du_theta^(-1)*k0*cN05^(-1)*kron((Isym*dNx(mapdNx,:))',N_G(k,:));
        kTc = N(k,:)'*(Deta_theta)^(-1)*D2eta_cc*(DotC(map.Voigt)'*(Galgo(map.Voigt)'*Isym)')*N_c(k,:) - k0*cN05^(-2)*(Du_theta)^(-1)*(dNx'*(GN05*(dNx*thetaN05')))*N_c(k,:);
        DQHN05_thetaN05 = -k0*cN05^(-1)*GN05*dNx*eye(size(dNx,2));
        kTT = 2*N(k,:)'*1/DT*N(k,:) - N(k,:)'*Deta_theta^(-2)*(D2eta_thetatheta*Deta_c)*N(k,:)*(sum(sum(DotC.*Galgo))) - Du_theta^(-1)*dNx'*DQHN05_thetaN05;
        KTEPS = KTEPS + 0.5*[2*kTC kTG kTc]*detJ*wp(k);
        KTT = KTT + 0.5*kTT*detJ*wp(k);

        ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 DeltaU];
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
