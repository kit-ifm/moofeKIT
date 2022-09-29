function out = dispCascadeSC_mooneyRivlin_midPoint(obj,varargin)
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

flagNumericalTangent = true;
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
    numberOfTemperatureDOFs=numel(temperatureDOFs);
    numberOfMechanicalDOFs=numel(mechanicalDOFs);
    Re=zeros(numDOFs,1);
    Ke=zeros(numDOFs);
    
    %% Residual and tangent
    [RX,RT,KXX,KTT,ePot] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
    if flagNumericalTangent
        %% Numerical Tangent
        epsilon=10^-7;
        KXXxx = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
            KXXxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KXTxx = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT);
            KXTxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KTXxx = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
            KTXxx(:,xx) = (RTxx-RT)/epsilon;
        end
        KTTxx = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT);
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

function [RX,RT,KXX,KTT,ePot] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT)
map.Voigt = [1 5 9 4 8 7]';
%% Shape functions
N = objNT.SHAPEF.N;
wp = objNT.SHAPEF.wp;
NGP = objNT.NGP;
DIM = objNT.DIM;
%% Material parameters
a = objNT.MAT.c1;
b = objNT.MAT.c2;
c = objNT.MAT.c;
d = 2*(a+2*b);
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
KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
% Initialize sub-residual vector
RX = zeros(numberOfMechanicalDOFs,1);
RT = zeros(numberOfTemperatureDOFs,1);

edN05 = 0.5*(edN1+edN);
thetaN05 = 0.5*(thetaN1+thetaN);

ePot = 0;

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
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN05 = 0.5*wedge(CxN05,CxN05);
    cxN = det(CxN);
    cxN1 = det(CxN1);
    cxN05 = det(CxN05);    
    
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
    
    % Strain energy function
    PsiIsoN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    PsiVolN1 = -d*log(sqrt(cxN1)) + c/2*(sqrt(cxN1)-1)^2;
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c*(sqrt(cxN1) - 1) - d/sqrt(cxN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    ePot = ePot + PsiN1*detJ*wp(k);
    
    %% Impuls balance
    % Derivative of the strain energy function
    DPsi_C = a*eye(3);
    DPsi_G = b*eye(3);
    DPsi_cxN05 = -d/(2*cxN05) + c/2*(1-1/sqrt(cxN1)) - DIM*beta*(thetaN1e - theta0)*(c/2*1/sqrt(cxN05) + d/2*(cxN05)^(-3/2));

    S_CN05 = 2*(DPsi_C + wedge(DPsi_G,CxN05) + DPsi_cxN05*GxN05);        
    rX = BN05'*0.5*S_CN05(map.Voigt);
    
    etaN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c*(sqrt(cxN1)-1)-d/sqrt(cxN1));
    etaN = kappa*log(thetaNe/theta0)+DIM*beta*(c*(sqrt(cxN)-1)-d/sqrt(cxN));
    QN05 = -k0*(1/cxN05*GxN05*(dNx*thetaN05'));
    rT = (N(k,:)'*thetaN05e*(etaN1-etaN)/DT-dNx'*QN05);

    %% summation of tangent and residual
    RX = RX + rX*detJ*wp(k);
    RT = RT + rT*detJ*wp(k);
        
end
end
