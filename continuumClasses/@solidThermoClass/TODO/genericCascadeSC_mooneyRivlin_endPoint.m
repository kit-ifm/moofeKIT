function out = genericCascadeSC_mooneyRivlin_endPoint(obj,varargin)
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
LocalDofs = [];
stress = 0;

for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'LocalDofs')
        LocalDofs = varargin{i+1};
    elseif strcmp(varargin{i},'stress')
        stress = varargin{i+1};
    end
end

dNrAll = obj.SHAPEF.dNr;
% Element stuff
edof = obj.EDOF;
numElements = numel(edof);
edof = obj.GLOBALEDOF;
edofZw = obj.EDOF;
QREF = obj.QREF;
QN1 = obj.QN1;
QN = obj.QN;
map.Voigt = [1 5 9 4 8 7]';
map.VoigtInv = [1 4 6; 4 2 5; 6 5 3];
map.VoigtFull = [1 5 9 4 8 7 2 6 3]';
map.VoigtFullInv = [1 4 6; 7 2 5; 9 8 3];
DIM = obj.DIM;
%% Shape functions
N = obj.SHAPEF.N;
wp = obj.SHAPEF.wp;
NGP = obj.NGP;
%% Material parameters
a = obj.MAT.c1;
b = obj.MAT.c2;
c1 = obj.MAT.c;
d1 = obj.MAT.d;
c2 = obj.MAT.cc;
d2 = obj.MAT.dd;
kappa = obj.MAT.kappa;
beta = obj.MAT.beta;
theta0 = obj.MAT.theta0;
k0 = obj.MAT.k0;
    
out(numElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
parfor j = 1:numElements
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
    % Number of external dofs
    numberOfTemperatureDOFs = numel(temperatureDOFs);
    numberOfMechanicalDOFs = numel(mechanicalDOFs);
    % Initialize sub-tangent matrices
%     KXX = zeros(numberOfMechanicalDOFs,numberOfMechanicalDOFs);
%     KXT = zeros(numberOfMechanicalDOFs,numberOfTemperatureDOFs);
%     KTX = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
%     KTT = zeros(numberOfTemperatureDOFs,numberOfTemperatureDOFs);
%     % Initialize sub-residual vector
%     RX = zeros(numberOfMechanicalDOFs,1);
%     RT = zeros(numberOfTemperatureDOFs,1);
    
    HelmholtzN1 = 0;
    EntropyN1 = 0;
    InternalEnergyN1 = 0;

    EntropyN = 0;
    EntropyDiff = 0;
    ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 EntropyDiff];
    
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
        
        % Deformation Gradient
        FxN = edN*dNx';
        FxN1 = edN1*dNx';
        CxN = FxN'*FxN;
        CxN1 = FxN1'*FxN1;
        DotCx = 1/DT*(CxN1-CxN);
        CxN = FxN'*FxN;
        GxN = 0.5*wedge(CxN,CxN);
        GxN1 = 0.5*wedge(CxN1,CxN1);
        cxN = det(CxN);
        cxN1 = det(CxN1);
        
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
        PsiVolN1 = c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1));
        PsiVolN = c1/2*(sqrt(cxN)-1)^2 - d1*log(sqrt(cxN));
        PsiThermoN1 = kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0));
        PsiCoupledN1 = -DIM*beta*(thetaN1e - theta0)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
        PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
        eThermN1 = kappa*(thetaN1e - theta0);
        eThermN = kappa*(thetaNe - theta0);
        eCoupledN1 = DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
        eCoupledN = DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
        entropyN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1)-d2/sqrt(cxN1));
        etotalN1 = PsiIsoN1 + PsiVolN1 + eThermN1 + eCoupledN1;

        entropyN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1)-d2/sqrt(cxN));
        
        HelmholtzN1 = HelmholtzN1 + PsiN1*detJ*wp(k);
        EntropyN1 = EntropyN1 + entropyN1*detJ*wp(k);
        InternalEnergyN1 = InternalEnergyN1 + etotalN1*detJ*wp(k);
        
        EntropyN = EntropyN + entropyN*detJ*wp(k);
        
        % MR-Mike
%         Psi = a*(trace(CN1)-DIM) + b*(1/2*((trace(CN1))^2 - trace((CN1)^2))-DIM) - 2*(a+2*b)*log((det(CN1))^0.5) + DIM*beta*(thetaN1 - theta0)*2*(a+2*b)*I3N1^(-1) + kappa*(thetaN1 - theta0 -thetaN1*log(thetaN1/theta0));
%         
%         PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - theta0 - thetaN1e*log(thetaN1e/theta0)) - DIM*beta*(thetaN1e - theta0)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
% Uebereinstimmung fuer:
%         c1 = 0;
%         c2 = 0;
%         d2 = d1;
        
        %% Impuls balance
        % Derivative of the strain energy function
        Dinternale_CN1 = a*eye(3);
        Dinternale_GN1 = b*eye(3);
        Dinternale_cxN1 = -d1/(2*cxN1) + c1/2*(1-1/sqrt(cxN1)) + DIM*beta*theta0*(c2/2*1/sqrt(cxN1) + d2/2*(cxN1)^(-3/2));
        Dentropy_cxN1 =  DIM*beta*(c2/2*1/sqrt(cxN1) + d2/2*(cxN1)^(-3/2));
        Dentropy_thetaN1 = kappa/thetaN1e;
        Dinternale_theta = kappa;
        D2internale_cxN1cxN1 = c1/(4*cxN1^(3/2)) + d1/(2*cxN1^2) - DIM*beta*theta0*(c2/(4*cxN1^(3/2)) + (3*d2)/(4*cxN1^(5/2)));
        D2entropy_cxN1cxN1 = -DIM*beta*(c2/(4*cxN1^(3/2)) + (3*d2)/(4*cxN1^(5/2)));
        
        % 2.P.-K. Stress tensor
        SN1 = 2*(Dinternale_CN1 + wedge(Dinternale_GN1,CxN1) + (Dinternale_cxN1-thetaN1e*Dentropy_cxN1)*GxN1);
        
        % Temperature evolution
        % Fourier law
        %QxN1 = -k0*(1/sqrt(cxN1)*GxN1*(dNx*thetaN1'));
        QxN1 = -k0*(1*cxN1^(-1)*GxN1*(dNx*thetaN1'));
        
        if stress~=0 % CALCULATE VON MISES STRESS FOR PLOT
            sigmaN1 = 1/det(FxN1)*FxN1*SN1*FxN1';
            tempStress = selectStress(sigmaN1,stress,DIM);
            
            Ke(temperatureDOFs,temperatureDOFs) = Ke(temperatureDOFs,temperatureDOFs) + (N(k,:)'*N(k,:))*detJ*wp(k);
            Re(temperatureDOFs,1) = Re(temperatureDOFs,1) + N(k,:)'*tempStress*detJ*wp(k);  
            
        else  % CALCULATE RESIDUAL AND TANGENT
            %% Residual densities
            rX = BN1'*0.5*SN1(map.Voigt);
            rT = N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT-(Dinternale_theta)^(-1)*dNx'*QxN1+N(k,:)'*(Dentropy_thetaN1)^(-1)*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*Isym)'*Dentropy_cxN1;

            %% Tangent operators densities
            % Lineraization of wedge(Dinternale_GN1,CxN1)
            D=(Dinternale_GN1);
            SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0  ;
                D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
                D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
                0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
                -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
                0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
            SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
            SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
            Kmat1 = SecDiffOperator;

            % Lineraization of (Dinternale_cxN1-thetaN1e*Dentropy_cxN1)*GxN1  part I
            D=(CxN1);
            SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0   ;
                D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
                D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
                0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
                -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
                0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
            SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
            SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
            Kmat2 = (Dinternale_cxN1-thetaN1e*Dentropy_cxN1)*SecDiffOperator;

            % Lineraization of (Dinternale_cxN1-thetaN1e*Dentropy_cxN1)*GxN1  part II
            Kmat3 = (D2internale_cxN1cxN1-thetaN1e*D2entropy_cxN1cxN1)*(GxN1(map.Voigt)*GxN1(map.Voigt)');

            % Material part of tangent operator
            ElasticityTensor = 2*Kmat1  + 2*Kmat2 + 4*Kmat3;
            TangentMaterial = 0.25*BN1'*ElasticityTensor*BN1;

            % Geometrical part of tangent operator
            GeometricalpartSave = dNx'*SN1*dNx;
            TangentGeometrical = zeros(numberOfMechanicalDOFs);
            for g = 1:DIM
                TangentGeometrical(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = GeometricalpartSave;
            end
            % Assembly of KXX
            kXX = TangentMaterial + TangentGeometrical;

            % KXTheta
            kXT = kron(BN1'*(-Dentropy_cxN1)*GxN1(map.Voigt),N(k,:));

            % KThetaThetaSecDiffOperator
            D2entropy_thetaN1thetaN1 = -kappa*thetaN1e^(-2)*N(k,:);
            %DQN1_thetaN1 = -k0*(1*cxN1^(-1/2)*GxN1*dNx*eye(size(dNx,2)));
            DQN1_thetaN1 = -k0*(1*cxN1^(-1)*GxN1*dNx*eye(size(dNx,2)));
            kTT = N(k,:)'*1/DT*N(k,:) - Dinternale_theta^(-1)*dNx'*DQN1_thetaN1 - N(k,:)'*Dentropy_thetaN1^(-2)*D2entropy_thetaN1thetaN1*(DotCx(map.Voigt)'*(GxN1(map.Voigt)'*Isym)'*Dentropy_cxN1);

            % KThetaX
            kTX1 = N(k,:)'*(Dentropy_thetaN1)^(-1)*DotCx(map.Voigt)'*(GxN1(map.Voigt)'*Isym)'*D2entropy_cxN1cxN1*(BN1'*GxN1(map.Voigt))';
            temp = wedge(CxN1,DotCx);
            kTX2 = N(k,:)'*(Dentropy_thetaN1)^(-1)*Dentropy_cxN1*(BN1'*temp(map.Voigt))';
            kTX3 = N(k,:)'*(Dentropy_thetaN1)^(-1)*Dentropy_cxN1*(1/DT*BN1'*GxN1(map.Voigt))';
            %kTX4 = -1/2*dNx'*Dinternale_theta^(-1)*k0*cxN1^(-3/2)*(GxN1*(dNx*thetaN1'))*(BN1'*GxN1(map.Voigt))';
            kTX4 = -dNx'*Dinternale_theta^(-1)*k0*cxN1^(-2)*(GxN1*(dNx*thetaN1'))*(BN1'*GxN1(map.Voigt))';
            kTX5 = zeros(numberOfTemperatureDOFs,numberOfMechanicalDOFs);
            for jj = 1:numberOfMechanicalDOFs
                BN1jj = Isyminv*BN1(:,jj);
                %kTX5(:,jj) = dNx'*Dinternale_theta^(-1)*k0*cxN1^(-1/2)*wedge(CxN1,BN1jj(map.VoigtInv))*(dNx*thetaN1');
                kTX5(:,jj) = dNx'*Dinternale_theta^(-1)*k0*cxN1^(-1)*wedge(CxN1,BN1jj(map.VoigtInv))*(dNx*thetaN1');
            end

            kTX = kTX1 + kTX2 + kTX3 + kTX4 + kTX5;

            %% summation of tangent and residual
%             RX = RX + rX*detJ*wp(k);
%             RT = RT + rT*detJ*wp(k);
% 
%             KXX = KXX + kXX*detJ*wp(k);
%             KXT = KXT + kXT*detJ*wp(k);
%             KTX = KTX + kTX*detJ*wp(k);
%             KTT = KTT + kTT*detJ*wp(k);
            Re(mechanicalDOFs,1) = Re(mechanicalDOFs,1) + rX*detJ*wp(k);
            Re(temperatureDOFs,1) = Re(temperatureDOFs,1) + rT*detJ*wp(k);
            
            Ke(mechanicalDOFs,mechanicalDOFs) = Ke(mechanicalDOFs,mechanicalDOFs) + kXX*detJ*wp(k);
            Ke(mechanicalDOFs,temperatureDOFs) = Ke(mechanicalDOFs,temperatureDOFs) + kXT*detJ*wp(k);
            Ke(temperatureDOFs,mechanicalDOFs) = Ke(temperatureDOFs,mechanicalDOFs) + kTX*detJ*wp(k);
            Ke(temperatureDOFs,temperatureDOFs) = Ke(temperatureDOFs,temperatureDOFs) + kTT*detJ*wp(k);
            
            ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 EntropyN1-EntropyN];
        end
    end
%     Re(mechanicalDOFs,1) = RX;
%     Re(temperatureDOFs,1) = RT;
%     Ke(mechanicalDOFs,mechanicalDOFs) = KXX;
%     Ke(mechanicalDOFs,temperatureDOFs) = KXT;
%     Ke(temperatureDOFs,mechanicalDOFs) = KTX;
%     Ke(temperatureDOFs,temperatureDOFs) = KTT;
    out(j).Re = Re;
    out(j).pK = Ke(:);
    out(j).ePot = ePot;
end
