function out = genericCascadeSC_mooneyRivlin_discreteGradientTCNumTang(obj,varargin)
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

flagNumericalTangent = true;
% Element stuff
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];


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
     flagDiscreteGradient = false(obj.NGP,2);
    [RX,RT,KXX,KTT,ePot,flagDiscreteGradient] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,0);
    if flagNumericalTangent
        %% Numerical Tangent
        epsilon=10^-7;
        KXXxx = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,flagNumericalTangent);
            KXXxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KXTxx = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [RXxx,~,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,flagNumericalTangent);
            KXTxx(:,xx) = (RXxx-RX)/epsilon;
        end
        KTXxx = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
        for xx = 1:1:numberOfMechanicalDOFs
            edN1xx = edN1;
            edN1xx(xx) = edN1xx(xx)+epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numDOFs,edN,edN1xx,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,flagNumericalTangent);
            KTXxx(:,xx) = (RTxx-RT)/epsilon;
        end
        KTTxx = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
        for xx = 1:1:numberOfTemperatureDOFs
            thetaN1xx = thetaN1;
            thetaN1xx(xx) = thetaN1xx(xx) + epsilon;
            [~,RTxx,~,~,~] = Residuum(DT,numDOFs,edN,edN1,thetaN1xx,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,flagNumericalTangent);
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

function [RX,RT,KXX,KTT,ePot,flagDiscreteGradient] = Residuum(DT,numDOFs,edN,edN1,thetaN1,thetaN,QREF,edofLocal,dNr,objNT,map,flagDiscreteGradient,flagNumericalTangent)
map.Voigt = [1 5 9 4 8 7]';
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
KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
% Initialize sub-residual vector
RX = zeros(numberOfMechanicalDOFs,1);
RT = zeros(numberOfTemperatureDOFs,1);

edN05 = 0.5*(edN1+edN);
thetaN05 = 0.5*(thetaN1+thetaN);

  HelmholtzN1 = 0;
    EntropyN1 = 0;
    InternalEnergyN1 = 0;
    DeltaInternale = 0;
    
ePot = [0 0 0 0];
Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
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
    DotCx = 1/DT*(CxN1-CxN);
    CxN05 = FxN05'*FxN05;
    CxAverageN05 = 0.5*(CxN1+CxN);
    GxN = 0.5*wedge(CxN,CxN);
    GxN1 = 0.5*wedge(CxN1,CxN1);
    GxN05 = 0.5*wedge(CxN05,CxN05);
    GxAverageN05 = 0.5*(GxN1+GxN);
    Unrat =1/3*(wedge(CxAverageN05,CxAverageN05) + GxAverageN05);
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
    PsiIsoN = a*(trace(CxN)-3) + b*(trace(GxN)-3);
    PsiVolN1 = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2;
    PsiVolN = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2;
    PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
    PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    eThermN1 = kappa*(thetaN1e - theta0);
    eThermN = kappa*(thetaNe - theta0);
    eCplN1 = DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    eCplN = DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
    etaN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
    etotalN1 = PsiIsoN1 + PsiVolN1 + eThermN1 + eCplN1;
    etotalN = PsiIsoN + PsiVolN + eThermN + eCplN;

    HelmholtzN1 = HelmholtzN1 + PsiN1*detJ*wp(k);
    EntropyN1 = EntropyN1 + etaN1*detJ*wp(k);
    InternalEnergyN1 = InternalEnergyN1 + etotalN1*detJ*wp(k);
    DeltaInternale = DeltaInternale + (etotalN1 - etotalN)*detJ*wp(k);
    %% Residual
    % Impuls evolution
    % Derivative of the strain energy function
    De_C = a*eye(3);
    De_G = b*eye(3);

    Dinternale_theta = kappa;
    % partitioned greenspan formulas
    if ~flagNumericalTangent
        
        if abs(thetaN1e-thetaNe) >= 10*eps
            flagDiscreteGradient(k,2) = true;
            % Entropy
            % N1N: mechanical field at N1 and theta at N
            entropyN1N1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyN1N = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyNN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            entropyNN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            Dentropy_thetaLeft = (entropyN1N1 -entropyN1N )*(thetaN1e-thetaNe)^(-1);
            Dentropy_thetaRight = (entropyNN1 -entropyNN )*(thetaN1e-thetaNe)^(-1);
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            theta =  Dentropy_theta^(-1)*Dinternale_theta;
        else
            % Entropy
            entropy_tempNN05 = kappa/thetaN05e;
            entropy_tempN1N05 = kappa/thetaN05e;
            Dentropy_thetaLeft = entropy_tempNN05;
            Dentropy_thetaRight = entropy_tempN1N05;
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            theta =  Dentropy_theta^(-1)*Dinternale_theta;
        end
        if abs(cxN1 - cxN) >= 10*eps
            flagDiscreteGradient(k,1) = true;
            % Energy
            % ec: all internal energy densities depending on c
            ecN1N1 = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2+kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecN1N = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2+kappa*(thetaNe - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecNN1 = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2+kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            ecNN = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2+kappa*(thetaNe - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            De_cxLeft=(ecN1N-ecNN)/(cxN1 - cxN);
            De_cxRight=(ecN1N1-ecNN1)/(cxN1 - cxN);
            De_cx= 0.5*(De_cxLeft+De_cxRight);
            % Entropy
            entropyN1N1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyN1N = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyNN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            entropyNN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            Dentropy_cxLeft=(entropyN1N-entropyNN)/(cxN1 - cxN);
            Dentropy_cxRight=(entropyN1N1-entropyNN1)/(cxN1 - cxN);
            Dentropy_cx= 0.5*(Dentropy_cxLeft+Dentropy_cxRight);
            % 2.P.K
            S = 2*(De_C + wedge(De_G,CxAverageN05) +  1/3*(De_cx-theta*Dentropy_cx)*(wedge(CxAverageN05,CxAverageN05) + GxAverageN05));
        else
            % Energy
            De_cx = -d1/(2*cxN05) + c1/2*(1-1/sqrt(cxN05)) + DIM*beta*theta0*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
            % Entropy
            Dentropy_cx =  DIM*beta*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
            % 2.P.K
            S = 2*(De_C + wedge(De_G,CxN05) + (De_cx-theta*Dentropy_cx)*GxN05);
        end
    else
        if  flagDiscreteGradient(k,2) == true
            % Entropy
            %N1N: mechanical field at N1 and theta at N
            entropyN1N1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyN1N = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyNN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            entropyNN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            Dentropy_thetaLeft = (entropyN1N1 -entropyN1N )*(thetaN1e-thetaNe)^(-1);
            Dentropy_thetaRight = (entropyNN1 -entropyNN )*(thetaN1e-thetaNe)^(-1);
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            theta =  Dentropy_theta^(-1)*Dinternale_theta;
        else
            % Entropy
            entropy_tempNN05 = kappa/thetaN05e;
            entropy_tempN1N05 = kappa/thetaN05e;
            Dentropy_thetaLeft = entropy_tempNN05;
            Dentropy_thetaRight = entropy_tempN1N05;
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            theta =  Dentropy_theta^(-1)*Dinternale_theta;
        end
        if flagDiscreteGradient(k,1) == true
            % Energy
            % ec: all internal energy densities depending on c
            ecN1N1 = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2+kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecN1N = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2+kappa*(thetaNe - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecNN1 = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2+kappa*(thetaN1e - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            ecNN = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2+kappa*(thetaNe - theta0) + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            De_cxLeft=(ecN1N-ecNN)/(cxN1 - cxN);
            De_cxRight=(ecN1N1-ecNN1)/(cxN1 - cxN);
            De_cx= 0.5*(De_cxLeft+De_cxRight);
            % Entropy
            entropyN1N1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyN1N = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
            entropyNN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            entropyNN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
            Dentropy_cxLeft=(entropyN1N-entropyNN)/(cxN1 - cxN);
            Dentropy_cxRight=(entropyN1N1-entropyNN1)/(cxN1 - cxN);
            Dentropy_cx= 0.5*(Dentropy_cxLeft+Dentropy_cxRight);
            % 2.P.K
            S = 2*(De_C + wedge(De_G,CxAverageN05) +  1/3*(De_cx-theta*Dentropy_cx)*(wedge(CxAverageN05,CxAverageN05) + GxAverageN05));
        else
            % Energy
            De_cx = -d1/(2*cxN05) + c1/2*(1-1/sqrt(cxN05)) + DIM*beta*theta0*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
            % Entropy
            Dentropy_cx =  DIM*beta*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
            % 2.P.K
            S = 2*(De_C + wedge(De_G,CxN05) + (De_cx-theta*Dentropy_cx)*GxN05);
        end
    end
    % Temperature evolution
    % Fourier law 
    QN05 = -k0*(cxN05^(-1)*GxAverageN05*(dNx*thetaN05'));
    %% Residual densities
    % impuls evolution
    rX = BN05'*0.5*S(map.Voigt);
    % temperature evolution
    rT =  (N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT-(Dinternale_theta)^(-1)*dNx'*QN05+N(k,:)'*(Dentropy_theta)^(-1)* DotCx(map.Voigt)'*(Unrat(map.Voigt)'*Isym)'*Dentropy_cx);
    %% summation of tangent and residual
    RX = RX + rX*detJ*wp(k);
    RT = RT + rT*detJ*wp(k);
        
    ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 DeltaInternale];
end
end
