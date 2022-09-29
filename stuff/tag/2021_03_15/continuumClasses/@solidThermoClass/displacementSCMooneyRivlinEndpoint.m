function out = displacementSCMooneyRivlinEndpoint(obj,setupObject)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = dispCascadeSC_mooneyRivlinThermo_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% dispCascadeSC=: displacement based cascPade formulation
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 05.04.2017 M.S., A.J., M.F.

DT = setupObject.timeStepSize;
globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
dimension = obj.dimension;
flagNumericalTangent = true;
objNT = obj;

out(numberOfElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
zw  = struct('RSe',[],'KSe',[]);

for e = 1:numberOfElements
    [edofH1,edofH2] = expandEdof(globalFullEdof(e,:));
    out(e).pI = edofH1';
    out(e).pJ = edofH2';
    out(e).edofE = double(globalFullEdof(e,:))';
    edN = qN(edof(e,:),1:dimension)';
    edN1 = qN1(edof(e,:),1:dimension)';
    thetaN = qN(edof(e,:),dimension+1)';
    thetaN1 = qN1(edof(e,:),dimension+1)';
    % External dofs
    temperatureDOFs = 4:dimension+1:numberOfDOFs;
    mechanicalDOFs = 1:numberOfDOFs;
    mechanicalDOFs(4:dimension+1:numberOfDOFs) = [];
    % Number of external dofs
    numberOfTemperatureDOFs=numel(temperatureDOFs);
    numberOfMechanicalDOFs=numel(mechanicalDOFs);
    Re=zeros(numberOfDOFs,1);
    Ke=zeros(numberOfDOFs);
    
    %% Residual and tangent
    [RX,RT,KXX,KTT,ePot] = Residuum(DT,numberOfDOFs,edN,edN1,thetaN1,thetaN,qR,edof(e,:),objNT);
    if flagNumericalTangent
        %% Numerical Tangent
        epsilon=10^-7;
        KXXxx = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numberOfDOFs,edN,edN1xx,thetaN1,thetaN,qR,edof(e,:),objNT);
            KXXxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KXTxx = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numberOfDOFs,edN,edN1,thetaN1xx,thetaN,qR,edof(e,:),objNT);
            KXTxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KTXxx = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numberOfDOFs,edN,edN1xx,thetaN1,thetaN,qR,edof(e,:),objNT);
            KTXxx(:,xx) = (RTxx-RT)/epsilon;
        end
        KTTxx = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numberOfDOFs,edN,edN1,thetaN1xx,thetaN,qR,edof(e,:),objNT);
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
    out(e).Re = Re;
    out(e).pK = Ke(:);
    out(e).ePot = ePot;
end
% obj.StaConStruct = zw;
end

function [RX,RT,KXX,KTT,ePot] = Residuum(DT,numberOfDOFs,edN,edN1,thetaN1,thetaN,qR,edofE,obj)
map.Voigt = [1 5 9 4 8 7]';
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
NAll = obj.shapeFunctions.N';
dNrAll = obj.shapeFunctions.dNr';
dimension = obj.dimension;

a = obj.materialData.a;
b = obj.materialData.b;
c = obj.materialData.c;
d = 2*(a+2*b);
kappa = obj.materialData.kappa;
beta = obj.materialData.beta;
theta0 = obj.materialData.theta0;
k0 = obj.materialData.k0;
% External dofs
temperatureDOFs = 4:dimension+1:numberOfDOFs;
mechanicalDOFs = 1:numberOfDOFs;
mechanicalDOFs(4:dimension+1:numberOfDOFs) = [];
% Number of external dofs
numberOfTemperatureDOFs = numel(temperatureDOFs);
numberOfMechanicalDOFs = numel(mechanicalDOFs);
% Initialize sub-tangent matrices
KXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
% Initialize sub-residual vector
RX = zeros(numberOfMechanicalDOFs,1);
RT = zeros(numberOfTemperatureDOFs,1);

ePot = 0;

J = qR(edofE,1:dimension)'*dNrAll;
% Run through all Gauss points
for k = 1:numberOfGausspoints
    indx = dimension*k-(dimension-1):dimension*k;
    detJ = det(J(:,indx)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,indx)')\(dNrAll(:,indx))';
%     dNx = (J(:,indx)')\dNr(indx,:);
    
    % Temperature
    thetaN1e = NAll(k,:)*thetaN1';
    thetaNe = NAll(k,:)*thetaN';
    
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
    CxN1 = FxN1'*FxN1;
    CxN = FxN'*FxN;
    GxN1 = 0.5*wedge(CxN1,CxN1);
    cxN1 = det(CxN1);    
    cxN = det(CxN);    
    
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
    PsiVolN1 = -d*log(sqrt(cxN1)) + c/2*(sqrt(cxN1)-1)^2;
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiCoupledN1 = -dimension*beta*(thetaN1e - theta0)*(c*(sqrt(cxN1) - 1) - d/sqrt(cxN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    ePot = ePot + PsiN1*detJ*gaussWeight(k);
    
    %% Impuls balance
    % Derivative of the strain energy function
    DPsi_C = a*eye(3);
    DPsi_G = b*eye(3);
    DPsi_cxN1 = -d/(2*cxN1) + c/2*(1-1/sqrt(cxN1)) - dimension*beta*(thetaN1e - theta0)*(c/2*1/sqrt(cxN1) + d/2*(cxN1)^(-3/2));
        
    S_CN1 = 2*(DPsi_C + wedge(DPsi_G,CxN1) + DPsi_cxN1*GxN1);        
    rX = BN1'*0.5*S_CN1(map.Voigt);
    
    etaN1 = kappa*log(thetaN1e/theta0)+dimension*beta*(c*(sqrt(cxN1)-1)-d/sqrt(cxN1));
    etaN = kappa*log(thetaNe/theta0)+dimension*beta*(c*(sqrt(cxN)-1)-d/sqrt(cxN));
    Q = -k0*(1/cxN1*GxN1*(dNx*thetaN1'));
    rT = (NAll(k,:)'*thetaN1e*(etaN1-etaN)/DT-dNx'*Q);

    %% summation of tangent and residual
    RX = RX + rX*detJ*gaussWeight(k);
    RT = RT + rT*detJ*gaussWeight(k);
        
end
end
