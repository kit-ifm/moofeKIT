function out = genericCascadeSC_mooneyRivlin_midPoint(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = dispCascadeSC_mooneyRivlinThermo_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% dispCascadeSC=: displacement based cascade formulation
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 05.04.2017 M.S., A.J., M.F.

%% Check input
DT = 1;
stress=false;
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'LocalDofs')
        LocalDofs = varargin{i+1};
    elseif strcmp(varargin{i},'stress')
        stress = true;
    end
end

flagNumericalTangent = false;
dNrAll = obj.SHAPEF.dNr;
objNT = obj;

% Element stuff
%% Create residual and tangent
edof = obj.EDOF;
numElements = numel(edof);
edof = obj.GLOBALEDOF;
edofZw = obj.EDOF;
QREF = obj.QREF;
QN1 = obj.QN1;
QN = obj.QN;
DIM = obj.DIM;
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];

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
    dNr = dNrAll;
    edN = QN(edofLocal,1:DIM)';
    edN1 = QN1(edofLocal,1:DIM)';
    thetaN = QN(edofLocal,DIM+1)';
    thetaN1 = QN1(edofLocal,DIM+1)';
    % External dofs
    temperatureDOFs = 4:DIM+1:numDOFs;
    mechanicalDOFs = 1:numDOFs;
    mechanicalDOFs(4:DIM+1:numDOFs) = [];
    % Number of external dofs
    Re=zeros(numDOFs,1);
    Ke=zeros(numDOFs);
    
    %% Residual and tangent
    [RX,RT,KXX,KXT,KTX,KTT,ePot] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,stress);
       if flagNumericalTangent
        %% Numerical Tangent
        epsilon=10^-7;
        numberOfTemperatureDOFs=numel(temperatureDOFs);
        numberOfMechanicalDOFs=numel(mechanicalDOFs);
        
        KXXxx = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [RXxx,~,~,~,~] =  Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,stress);
            KXXxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KXTxx = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [RXxx,~,~,~,~] =  Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT,map,stress);
            KXTxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KTXxx = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [~,RTxx,~,~,~] =  Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,stress);
            KTXxx(:,xx) = (RTxx-RT)/epsilon;
        end
         KTTxx = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT,map,stress);
            KTTxx(:,xx) = (RTxx-RT)/epsilon;
        end
        KXX = KXXxx;
        KXT = KXTxx;
        KTX = KTXxx;
         KTT = KTTxx;
    end
    Re(mechanicalDOFs,1) = RX;
    Re(temperatureDOFs,1) = RT;
    Ke(mechanicalDOFs,mechanicalDOFs) = KXX;
    Ke(mechanicalDOFs,temperatureDOFs) = KXT;
    Ke(temperatureDOFs,mechanicalDOFs) = KTX;
    Ke(temperatureDOFs,temperatureDOFs) = KTT;
    out(j).Re = Re;
    out(j).pK = Ke(:);
    out(j).ePot = ePot;
end
obj.StaConStruct = zw;
end

function [RX,RT,KXX,KXT,KTX,KTT,ePot] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,stress)
%% Shape functions
N = objNT.SHAPEF.N;
wp = objNT.SHAPEF.wp;
NGP = objNT.NGP;
DIM = objNT.DIM;
%% Material parameters
a = objNT.MAT.c1;
b = objNT.MAT.c2;
c1 = objNT.MAT.c;
d1 = objNT.MAT.d;
c2 = objNT.MAT.cc;
d2 = objNT.MAT.dd;
kappa = objNT.MAT.kappa;
beta = objNT.MAT.beta;
theta0 = objNT.MAT.theta0;
k0 = objNT.MAT.k0;
% External dofs
temperatureDOFs = 4:DIM+1:numDOFs;
mechanicalDOFs = 1:numDOFs;
mechanicalDOFs(4:DIM+1:numDOFs) = [];
% Number of external dofs
numberOfTemperatureDOFs = numel(temperatureDOFs);
numberOfMechanicalDOFs = numel(mechanicalDOFs);
% Initialize sub-tangent matrices
KXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
KXT = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
KTX = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
% Initialize sub-residual vector
RX = zeros(numberOfMechanicalDOFs,1);
RT = zeros(numberOfTemperatureDOFs,1);

HelmholtzN1 = 0;
EntropyN1 = 0;
InternalEnergyN1 = 0;
ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1];

edN05 = 0.5*(edN1+edN);
thetaN05 = 0.5*(thetaN1+thetaN);

I = eye(3);
Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
Isyminv = [eye(3) zeros(3); zeros(3) 1/2*eye(3)];
J = QREF(edofLocal,1:DIM)'*dNr';
% Run through all Gauss points
for k = 1:NGP
    indx = DIM*k-(DIM-1):DIM*k;
    detJ = det(J(:,indx)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,indx)')\dNr(indx,:);

    % Temperature
    thetaNe = N(k,:)*thetaN';
    thetaN1e = N(k,:)*thetaN1';
    thetaN05e = N(k,:)*thetaN05';
    
    % Deformation Gradient
    FxN = edN*dNx';
    FxN1 = edN1*dNx';
    FxN05 = edN05*dNx';
    CxN = FxN'*FxN;
    CxN1 = FxN1'*FxN1;
    CxN05 = FxN05'*FxN05;
    GxN = 0.5*wedge(CxN,CxN);
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN05 = 0.5*wedge(CxN05,CxN05);
    cxN = det(CxN);
    cxN1 = det(CxN1);
    cxN05 = det(CxN05);    
    DotCx = 1/DT*(CxN1-CxN);
    
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
    
      % Strain energy function
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN1 = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2;
    PsiVolN = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2;
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    emixN1 = kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    emixN = kappa*(thetaNe - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
    entropyN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
    entropyN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
    etotalN1 = PsiIsoN1 + PsiVolN1 + emixN1;
    etotalN = PsiIsoN + PsiVolN + emixN;
       
    HelmholtzN1 = HelmholtzN1 + PsiN1*detJ*wp(k);
    EntropyN1 = EntropyN1 + entropyN1*detJ*wp(k);
    InternalEnergyN1 = InternalEnergyN1 + etotalN1*detJ*wp(k);
    
    %% Impuls balance
    % Derivative of the strain energy function
    Dinternale_CN05 = a*eye(3);
    Dinternale_GN05 = b*eye(3);
    Dinternale_cxN05 = -d1/(2*cxN05) + c1/2*(1-1/sqrt(cxN05)) + DIM*beta*theta0*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));    
    Dentropy_cxN05 =  DIM*beta*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
    Dentropy_thetaN05 = kappa/thetaN05e;
    Dinternale_theta = kappa;
    D2internale_cxN05cxN05 = c1/(4*cxN05^(3/2)) + d1/(2*cxN05^2) - DIM*beta*theta0*(c2/(4*cxN05^(3/2))+ (3*d2)/(4*cxN05^(5/2)));
    D2entropy_cxcxN05 =-DIM*beta*(c2/(4*cxN05^(3/2)) + (3*d2)/(4*cxN05^(5/2)));
    
    % 2.P.-K. Stress tensor
    SN05 = 2*(Dinternale_CN05 + wedge(Dinternale_GN05,CxN05) + (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05);

    % Temperature evolution
    % Fourier law    
    QN05 = -k0*(1/cxN05*GxN05*(dNx*thetaN05'));
     if stress==true
            sigmaalgo = 1/det(FxN05)*FxN05*SN05*FxN05';
            tempStress = selectStress(sigmaalgo,-1,DIM);
            KTT =KTT + (N(k,:)'*N(k,:))*detJ*wp(k);
            RT = RT + N(k,:)'*tempStress*detJ*wp(k);  
        else    
    %% Residual densities
    rX = BN05'*0.5*SN05(map.Voigt);
    rT = N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT-(Dinternale_theta)^(-1)*dNx'*QN05+N(k,:)'*(Dentropy_thetaN05)^(-1)*DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*Dentropy_cxN05;

   %% Tangent operators densities
   % Lineraization of wedge(Dinternale_GN05,CxN05)    
    D=(Dinternale_GN05);
    SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0  ;
        D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
        D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
        0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
        -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
        0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
    SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
    SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
    Kmat1 = SecDiffOperator;
    
    % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part I
    D=(CxN05);
    SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0   ;
        D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
        D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
        0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
        -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
        0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
    SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
    SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
    Kmat2 = (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*SecDiffOperator;
    
    % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part II 
    Kmat3 = (D2internale_cxN05cxN05-thetaN05e*D2entropy_cxcxN05)*(GxN05(map.Voigt)*GxN05(map.Voigt)');    
    
    % Material part of tangent operator
    ElasticityTensor = 2*Kmat1  + 2*Kmat2 + 4*Kmat3;
    TangentMaterial = 0.25*BN05'*ElasticityTensor*BN05;
    
    % Geometrical part of tangent operator
    GeometricalpartSave = dNx'*SN05*dNx;
    TangentGeometrical = zeros(numberOfMechanicalDOFs);
    for g = 1:DIM
        TangentGeometrical(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = GeometricalpartSave;
    end
    % Assembly of KXX
    kXX = 0.5*(TangentMaterial + TangentGeometrical);
      
    % KXTheta
    kXT = 0.5*kron(BN05'*(-Dentropy_cxN05)*GxN05(map.Voigt),N(k,:));
    
    % KThetaThetaSecDiffOperator
    D2entropy_thetaN05thetaN05 = -kappa*thetaN05e^(-2)*N(k,:);
    DQN05_thetaN05 = -k0*(1/cxN05*GxN05*dNx*eye(size(dNx,2)));
    kTT = N(k,:)'*1/DT*N(k,:) - 0.5*Dinternale_theta^(-1)*dNx'*DQN05_thetaN05 - 0.5*N(k,:)'*Dentropy_thetaN05^(-2)*D2entropy_thetaN05thetaN05*(DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*Dentropy_cxN05);
    
    % KThetaX
    kTX1 = N(k,:)'*(Dentropy_thetaN05)^(-1)*DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*D2entropy_cxcxN05*(BN05'*GxN05(map.Voigt))';
    temp = wedge(CxN05,DotCx);
    kTX2 = N(k,:)'*(Dentropy_thetaN05)^(-1)*Dentropy_cxN05*(BN05'*temp(map.Voigt))';
    kTX3 = N(k,:)'*(Dentropy_thetaN05)^(-1)*Dentropy_cxN05*(1/DT*BN1'*GxN05(map.Voigt))';
    kTX4 = -dNx'*Dinternale_theta^(-1)*k0*cxN05^(-2)*(GxN05*(dNx*thetaN05'))*(BN05'*GxN05(map.Voigt))';    
    Q_1G=zeros(3,3);
    Q_2G=zeros(3,3);
    Q_3G=zeros(3,3);
    Q_1G(:,:)=-0.5*k0*1/cxN05*(I(:,1)*(dNx*thetaN05')'+(dNx*thetaN05')*I(1,:));
    Q_2G(:,:)=-0.5*k0*1/cxN05*(I(:,2)*(dNx*thetaN05')'+(dNx*thetaN05')*I(2,:));
    Q_3G(:,:)=-0.5*k0*1/cxN05*(I(:,3)*(dNx*thetaN05')'+(dNx*thetaN05')*I(3,:));
    Q_Gv = zeros(DIM,6);
    Q_Gv(1,:)=Q_1G(map.Voigt);
    Q_Gv(2,:)=Q_2G(map.Voigt);
    Q_Gv(3,:)=Q_3G(map.Voigt);
    D=(CxN05);
    SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0   ;
        D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
        D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
        0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
        -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
        0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
    SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
    SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
    kTX5 = -(Dinternale_theta)^(-1)*dNx'*Q_Gv*Isym*0.5*SecDiffOperator*BN05;
    kTX = 0.5*kTX1 + 0.5*kTX2 + kTX3 + 0.5*kTX4 + 0.5*kTX5;

    %% summation of tangent and residual
    RX = RX + rX*detJ*wp(k);
    RT = RT + rT*detJ*wp(k);
        
    KXX = KXX + kXX*detJ*wp(k);
    KXT = KXT + kXT*detJ*wp(k);
    KTX = KTX + kTX*detJ*wp(k);
    KTT = KTT + kTT*detJ*wp(k);
    
    ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1];
     end
end
end
