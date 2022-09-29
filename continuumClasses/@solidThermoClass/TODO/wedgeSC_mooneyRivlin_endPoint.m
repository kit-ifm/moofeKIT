function out = wedgeSC_mooneyRivlin_endPoint(obj,varargin)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = mmhwCascadeSC_mooneyRivlinThermo_discreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% wedgeSC=: displacment based formulation, formulated in a polyconvex
% framweork. Based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 05.04.2017 M.S., A.J., M.F.

% Check input
DT = 1;
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    end
end
ansatzfunction = obj.Ansatzfunction;
dNrAll = obj.SHAPEF.dNr;
numTang = true;
objNT = obj;
DIM = obj.DIM;

%% Create residual and tangent
edof = obj.GLOBALEDOF;
numElements = numel(edof);
edofZw = obj.EDOF;
QREF = obj.QREF;
QN1 = obj.QN1;
QN = obj.QN;

%% Parfor
out(numElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
for j = 1:numElements
    dNr = [];
    [edofH1,edofH2] = expandEdof(edof{j});
    out(j).pI = edofH1';
    out(j).pJ = edofH2';
    out(j).edofE = edof{j}';
    numDOFs = numel(edof{j});
    edofLocal = edofZw{j};
    if strcmpi(ansatzfunction,'Lagrange')
        dNr = dNrAll;
    elseif strcmpi(ansatzfunction,'HNURBS')
        dNr = dNrAll{j};
    end
    edN = QN(edofLocal,1:DIM)';
    edN1 = QN1(edofLocal,1:DIM)';
    thetaN1 = QN1(edofLocal,DIM+1)';
    thetaN = QN(edofLocal,DIM+1)';

    
    %% Residual and tangent
    [RX,RT,KXX,KXT,KTX,KTT,helmholtz,entropyN1,internaleN1] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
        
    % External dofs
        dofTemp=4:DIM+1:numDOFs;
        dofMech=1:numDOFs;
        dofMech(4:DIM+1:numDOFs)=[];
        
        % Number of external dofs
        numTemp=numel(dofTemp);
        numMech=numel(dofMech);
        
    if numTang == 1
    
        %% Numerical Tangent
        epsilon=10^-7;
        KXXxx = zeros(numMech,numMech);
        for xx = 1:1:numMech
            edN1_xx=edN1;
            edN1_xx(xx)=edN1_xx(xx)+epsilon;
            [RX_num,~,~,~,~,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1_xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
            KXXxx(:,xx) = (RX_num-RX)/epsilon;
        end
        KXTxx = zeros(numMech,numTemp);
        for xx = 1:1:numTemp
            thetaN1_xx=thetaN1;
            thetaN1_xx(xx)=thetaN1_xx(xx)+epsilon;
            [RX_num,~,~,~,~,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1_xx,thetaN,QREF,edofLocal,dNr,objNT);
            KXXxx(:,xx) = (RX_num-RX)/epsilon;
        end
        
        KTXtt = zeros(numTemp,numMech);
        for xx = 1:1:numMech
            edN1_xx=edN1;
            edN1_xx(xx)=edN1_xx(xx)+epsilon;
            [~,RT_num,~,~,~,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1_xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT);
            KTXtt(:,xx) = (RT_num-RT)/epsilon;
        end
        
        KTTxx = zeros(numTemp,numTemp);
        for xx = 1:1:numTemp
            thetaN1_xx=thetaN1;
            thetaN1_xx(xx)=thetaN1_xx(xx)+epsilon;
            [~,RT_num,~,~,~,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1_xx,thetaN,QREF,edofLocal,dNr,objNT);
            KTTxx(:,xx) = (RT_num-RT)/epsilon;
        end
    
        % Overwrite
        KXX = KXXxx;
        KXT = KXTxx;
        KTX = KTXtt;
        KTT = KTTxx;
    end
        Re(dofMech,1)=RX;
        Re(dofTemp,1)=RT;
        Ke(dofMech,dofMech)=KXX;
        Ke(dofMech,dofTemp)=KXT;
        Ke(dofTemp,dofMech)=KTX;
        Ke(dofTemp,dofTemp)=KTT;
        
        out(j).Re = Re;
        out(j).pK = Ke(:);
        out(j).ePot = [helmholtz entropyN1 internaleN1 ];
    end
end 


function [RX,RT,KXX,KXT,KTX,KTT,helmholtz,entropyN1,internaleN1]=Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT)

% Shape functions
wp = objNT.SHAPEF.wp;
NGP = objNT.NGP;
DIM = objNT.DIM;
N = objNT.SHAPEF.N;

% Material parameters
a = objNT.MAT.c1;
b = objNT.MAT.c2;
c = objNT.MAT.c;
if objNT.MAT.strongcomp == true
    d = 2*(a+2*b);
else
    d = 0;
end
kappa = objNT.MAT.kappa;
beta = objNT.MAT.beta;
theta0 = objNT.MAT.theta0;
k0 = objNT.MAT.k0;

% External dofs
dofTemp=4:DIM+1:numDOFs;
dofMech=1:numDOFs;
dofMech(4:DIM+1:numDOFs)=[];

% Number of external dofs
numTemp=numel(dofTemp);
numMech=numel(dofMech);

%Initialize sub-tangent matrices
KXX = zeros(numMech,numMech);
KXT = zeros(numMech,numTemp);
KTX = zeros(numTemp,numMech);
KTT = zeros(numTemp);
%Initialize sub-residual vector
RT = zeros(numTemp,1);
RX = zeros(numMech,1);
%Initialize  scalars
helmholtz = 0;
entropyN1 = 0;
internaleN1 = 0;

%Jacobi-Matrix
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
    thetaN1e = N(k,:)*thetaN1';
    thetaNe = N(k,:)*thetaN';
    
    % Deformation gradient, Right Cauchy-Green tensor
    FN1 = edN1*dNx';
    FN = edN*dNx';
    CN1 = FN1'*FN1;
    CN = FN'*FN;
    
    % B-matrix (current configuration)
    BN1 = zeros(6,numMech);
    BN1(1,1:3:end)=FN1(1,1)*dNx(1,:);
    BN1(1,2:3:end)=FN1(2,1)*dNx(1,:);
    BN1(1,3:3:end)=FN1(3,1)*dNx(1,:);
    BN1(2,1:3:end)=FN1(1,2)*dNx(2,:);
    BN1(2,2:3:end)=FN1(2,2)*dNx(2,:);
    BN1(2,3:3:end)=FN1(3,2)*dNx(2,:);
    BN1(3,1:3:end)=FN1(1,3)*dNx(3,:);
    BN1(3,2:3:end)=FN1(2,3)*dNx(3,:);
    BN1(3,3:3:end)=FN1(3,3)*dNx(3,:);
    BN1(4,1:3:end)=FN1(1,1)*dNx(2,:)+FN1(1,2)*dNx(1,:);
    BN1(4,2:3:end)=FN1(2,1)*dNx(2,:)+FN1(2,2)*dNx(1,:);
    BN1(4,3:3:end)=FN1(3,1)*dNx(2,:)+FN1(3,2)*dNx(1,:);
    BN1(5,1:3:end)=FN1(1,2)*dNx(3,:)+FN1(1,3)*dNx(2,:);
    BN1(5,2:3:end)=FN1(2,2)*dNx(3,:)+FN1(2,3)*dNx(2,:);
    BN1(5,3:3:end)=FN1(3,2)*dNx(3,:)+FN1(3,3)*dNx(2,:);
    BN1(6,1:3:end)=FN1(1,1)*dNx(3,:)+FN1(1,3)*dNx(1,:);
    BN1(6,2:3:end)=FN1(2,1)*dNx(3,:)+FN1(2,3)*dNx(1,:);
    BN1(6,3:3:end)=FN1(3,1)*dNx(3,:)+FN1(3,3)*dNx(1,:);
    BN1 = BN1;
    
    % Cofactor
    GN1 = 0.5*wedge(CN1,CN1);
    % Third invariant
    cN1 = det(CN1);
    cN = det(CN);
    
    % Strain energy function
    PsiIsoN1 = a*(trace(CN1)-3) + b*(trace(GN1)-3);
    PsiVolN1 = -d*log(sqrt(cN1)) + c/2*(sqrt(cN1)-1)^2;
    emechN1 = PsiVolN1;
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c*(sqrt(cN1) - 1) - d/sqrt(cN1));
    emixN1 = kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c*(sqrt(cN1) - 1) - d/sqrt(cN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    etotalN1 = PsiIsoN1 + emechN1 + emixN1;
    etaN1 = kappa*log(thetaN1e/theta0) + DIM*beta*(c*(sqrt(cN1)-1)-d/sqrt(cN1));
    helmholtz = helmholtz + PsiN1*detJ*wp(k);
    entropyN1 = entropyN1 + etaN1*detJ*wp(k);
    internaleN1 = internaleN1 + etotalN1*detJ*wp(k);
    
    % Impuls balance
    % Derivative of the strain energy function
    DPsi_C = a*eye(3);
    DPsi_G = b*eye(3);
    DPsi_cN1 = -d/(2*cN1) + c/2*(1-1/sqrt(cN1)) - DIM*beta*(thetaN1e - theta0)*(c/2*1/sqrt(cN1) + d/2*(cN1)^(-3/2));
    
    % Hessian operator of the complementary energy function
%     D2Psi_cc = d/(2*cN1^2) + c/(4*(cN1)^(3/2)) - DIM*beta*(thetaN1e - theta0)*(-c/4*cN1^(-3/2) - 3*d/4*cN1^(-5/2));
%     D2Psi_ct = DIM*beta*(c/2*1/sqrt(cN1)+d/2*(cN1)^(-3/2));
    
    % Second Piola Kirchhoff stress tensor
    SAkt = 2*(DPsi_C + wedge(DPsi_G,CN1) + 1/3*DPsi_cN1*(wedge(CN1,CN1)+GN1));
    SAkt_v = [SAkt(1,1); SAkt(2,2); SAkt(3,3); SAkt(1,2); SAkt(2,3); SAkt(1,3)];
    
    %% Energy Balance
    % Entropy
    etaN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c*(sqrt(cN1)-1)-d/sqrt(cN1));
    etaN = kappa*log(thetaNe/theta0)+DIM*beta*(c*(sqrt(cN)-1)-d/sqrt(cN));
    % Fourier law
    Q = -k0*(1/cN1*GN1*(dNx*thetaN1'));
    
    %% Residual
    RX = RX + BN1'*0.5*SAkt_v*detJ*wp(k);
    RT = RT + (N(k,:)'*thetaN1e*(etaN1-etaN)/DT-dNx'*Q)*detJ*wp(k);
    
%     %% Tangent
%     % Derivative of Sigma_I
%     Sigma_I_I = 2*d/I3N1^2+c/(I3N1)^(3/2); % S_I_I = 4*D2W/DI3^2
%     
%     % Derivative of GAkt
%     zwF = zeros(6);
%     zwF(2,1) = CN1(3,3);
%     zwF(3,1) = CN1(2,2);
%     zwF(3,2) = CN1(1,1);
%     zwF(4,3) = -CN1(2,1);
%     zwF(5,1) = -CN1(3,2);
%     zwF(6,4) = 0.5*CN1(3,2);
%     zwF(5,4) = 0.5*CN1(3,1);
%     zwF(6,2) = -CN1(3,1);
%     zwF(6,5) = 0.5*CN1(2,1);
%     zwF = 2*(zwF + zwF');
%     zwF(4,4) = -CN1(3,3);
%     zwF(5,5) = -CN1(1,1);
%     zwF(6,6) = -CN1(2,2);
%     
%     % Derivative of wedge(Sigma_G,CAkt)
%     zwF2 = zeros(6);
%     zwF2(2,1) = 2*c2;
%     zwF2(3,1) = 2*c2;
%     zwF2(3,2) = 2*c2;
%     zwF2 = 2*(zwF2 + zwF2');
%     zwF2(4,4) = -2*c2;
%     zwF2(5,5) = -2*c2;
%     zwF2(6,6) = -2*c2;
%     
%     % Assembly of elasticity tensor
%     GAkt_v = [GN1(1,1); GN1(2,2); GN1(3,3); GN1(1,2); GN1(2,3); GN1(1,3)];
%     ELA = zwF2 + Sigma_I*zwF + Sigma_I_I*(GAkt_v*GAkt_v');
%     
%     A1 = dNx'*SAkt*dNx*detJ*wp(k);
%     MAT = zeros(numDOFs);
%     for g = 1:DIM
%         MAT(g:DIM:numDOFs,g:DIM:numDOFs) = A1;
%     end
%     Ke = Ke + 0.25*BN1'*ELA*BN1*detJ*wp(k) + MAT;
end
end
