function out = genericCascadeSC_mooneyRivlin_discreteGradient(obj,varargin)
%% Creates the residual and the tangent of the given obj.
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
stress = false;
for i = 1:size(varargin,2)
    if strcmp(varargin{i},'dt')
        if ~isnumeric(varargin{i+1})
            error('Please provide a numeric as input.')
        end
        DT = varargin{i+1};
    elseif strcmp(varargin{i},'LocalDofs')
        LocalDofs = varargin{i+1};
    elseif  strcmp(varargin{i},'stress')
        stress=true;
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

Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];
Isyminv = [eye(3) zeros(3); zeros(3) 1/2*eye(3)];

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
    edN05 = 0.5*(edN1+edN);
    thetaN05 = 0.5*(thetaN1+thetaN);
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
    DeltaInternale = 0;
    ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 DeltaInternale];
    
    flagDiscreteGradient = false(obj.NGP,2);
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
        eCplN1 = DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
        eCplN = DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
        etaN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1) - d2/sqrt(cxN1));
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
        
        De_theta = kappa;
        Dentropy_cxN05 =  DIM*beta*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2));
        D2entropy_cxcxN05 = -DIM*beta*(c2/(4*cxN05^(3/2)) + (3*d2)/(4*cxN05^(5/2)));
        
        % partitioned greenspan formulas
        
        if abs(thetaN1e-thetaNe) >= 10*eps
            flagDiscreteGradient(k,2) = true;
            % Entropy
            entropyN1N1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1) - d2/sqrt(cxN1));
            entropyN1N = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN1)-1) - d2/sqrt(cxN1));
            entropyNN1 = kappa*log(thetaN1e/theta0)+DIM*beta*(c2*(sqrt(cxN)-1) - d2/sqrt(cxN));
            entropyNN = kappa*log(thetaNe/theta0)+DIM*beta*(c2*(sqrt(cxN)-1) - d2/sqrt(cxN));
            Dentropy_thetaLeft = (entropyN1N1 -entropyN1N )*(thetaN1e-thetaNe)^(-1);
            Dentropy_thetaRight = (entropyNN1 -entropyNN )*(thetaN1e-thetaNe)^(-1);
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            thetaAlgo =  Dentropy_theta^(-1)*De_theta;
            DentropyN1N1_theta =kappa/thetaN1e;
            DentropyNN1_theta = kappa/thetaN1e;
            D2entropy_thetatheta = 1/2*((DentropyN1N1_theta + DentropyNN1_theta)/(thetaN1e - thetaNe) - (thetaN1e - thetaNe)^(-2)*(entropyN1N1 - entropyN1N + entropyNN1 - entropyNN));
            %D2entropy_thetatheta2 = 2*D2entropy_thetatheta;
        else
            % Entropy
            entropy_tempNN05 = kappa/thetaN05e;
            entropy_tempN1N05 = kappa/thetaN05e;
            Dentropy_thetaLeft = entropy_tempNN05;
            Dentropy_thetaRight = entropy_tempN1N05;
            Dentropy_theta = 1/2*(Dentropy_thetaLeft+Dentropy_thetaRight);
            % Temperature
            thetaAlgo =  Dentropy_theta^(-1)*De_theta;
            D2entropy_thetatheta = -kappa*thetaN05e^(-2);
        end
        if abs(cxN1 - cxN) >= 10*eps
            flagDiscreteGradient(k,1) = true;
            % Energy
            % ec: all internal energy densities depending on c
            ecN1N1 = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2 + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecN1N = -d1*log(sqrt(cxN1)) + c1/2*(sqrt(cxN1)-1)^2 + DIM*beta*theta0*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
            ecNN1 = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2 + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            ecNN = -d1*log(sqrt(cxN)) + c1/2*(sqrt(cxN)-1)^2 + DIM*beta*theta0*(c2*(sqrt(cxN) - 1) - d2/sqrt(cxN));
            DecN1N_cx = - d1/(2*cxN1) + (c1*(cxN1^(1/2) - 1))/(2*cxN1^(1/2)) + DIM*beta*theta0*(c2/(2*cxN1^(1/2)) + d2/2*(cxN1)^(-3/2));
            DecN1N1_cx = - d1/(2*cxN1) + (c1*(cxN1^(1/2) - 1))/(2*cxN1^(1/2)) + DIM*beta*theta0*(c2/(2*cxN1^(1/2)) + d2/2*(cxN1)^(-3/2));
            Dec_cxLeft=(ecN1N-ecNN)/(cxN1 - cxN);
            Dec_cxRight=(ecN1N1-ecNN1)/(cxN1 - cxN);
            Dec_cx= 0.5*(Dec_cxLeft+Dec_cxRight);
            D2ec_DcxDcx= ((DecN1N_cx + DecN1N1_cx)/(cxN1 - cxN) + (-(cxN1 - cxN)^(-2))*(ecN1N1 - ecNN1 + ecN1N - ecNN));
            SAlgo = 2*(De_C + wedge(De_G,CxAverageN05) +  Dec_cx*Unrat-thetaAlgo*Dentropy_cxN05*GxN05);
        else
            %Energy
            Dec_cx = -d1/(2*cxN05) + c1/2*(1-1/sqrt(cxN05)) + DIM*beta*theta0*(c2/2*1/sqrt(cxN05) + d2/2*(cxN05)^(-3/2)); 
            D2ec_DcxDcx = c1/(4*cxN05^(3/2)) + d1/(2*cxN05^2) - DIM*beta*theta0*(c2/(4*cxN05^(3/2))+ (3*d2)/(4*cxN05^(5/2)));            
            % 2.P.K
            SAlgo = 2*(De_C + wedge(De_G,CxN05) + (Dec_cx-thetaAlgo*Dentropy_cxN05)*GxN05);
        end
        % Temperature evolution
        % Fourier law
        if stress==true
            sigmaalgo = 1/det(FxN05)*FxN05*SAlgo*FxN05';
            tempStress = selectStress(sigmaalgo,-1,DIM);
            KTT =KTT + (N(k,:)'*N(k,:))*detJ*wp(k);
            RT = RT + N(k,:)'*tempStress*detJ*wp(k);  
        else
        QN05 = -k0*(cxN05^(-1)*GxN05*(dNx*thetaN05'));
        %% Residual densities
        % impuls evolution
        rX = BN05'*0.5*SAlgo(map.Voigt);
        % temperature evolution
        rT =  (N(k,:)'*N(k,:)*(thetaN1-thetaN)'/DT-(De_theta)^(-1)*dNx'*QN05+N(k,:)'*(Dentropy_theta)^(-1)*DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*Dentropy_cxN05);
        %% summation residual
        RX = RX + rX*detJ*wp(k);
        RT = RT + rT*detJ*wp(k);
        %% Tangent
        % Lineraization of wedge(Dinternale_GN05,CxN05)
        D=(De_G);
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0  ;
            D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
            D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
            0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
            -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
            0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
        SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
        SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
        Kmat1 = SecDiffOperator;
        
        %    Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part I
        if flagDiscreteGradient(k,1) == false
            D=(CxN05);
            SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0   ;
                D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
                D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
                0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
                -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
                0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
            SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
            SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
            Kmat2 = Dec_cx*SecDiffOperator;
            
            % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part II
            Kmat3 = (D2ec_DcxDcx)*(GxN05(map.Voigt)*GxN05(map.Voigt)');
            
            % Material part of tangent operator
            ElasticityTensor = 2*Kmat1  + 2*Kmat2 + 4*Kmat3;
            TangentMateriale = 0.25*BN05'*ElasticityTensor*BN05;
        else
            % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part I
            D=(2/3*(CxAverageN05+0.5*CxN1));
            SecDiffOperator=[0          D(3,3)  D(2,2)      0 -D(3,2)  0   ;
                D(3,3)     0       D(1,1)      0         0       -D(3,1) ;
                D(2,2)     D(1,1)  0           -D(2,1)      0 0        ;
                0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3) ;
                -D(3,2)        0       0           D(3,1)   -D(1,1) D(1,2)   ;
                0          -D(3,1)    0           D(3,2)    D(2,1) -D(2,2)] ;
            SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
            SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
            Kmat2 = Dec_cx*SecDiffOperator;
            
            % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part II
            Kmat3 = (D2ec_DcxDcx)*Unrat(map.Voigt)*(GxN1(map.Voigt))';
            
            % Material part of tangent operator
            ElasticityTensor = 2*Kmat1  + 2*Kmat2 + 4*Kmat3;
            TangentMateriale = 0.25*BN05'*ElasticityTensor*BN1;
        end
        
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
        Kmat2 = (-thetaAlgo*Dentropy_cxN05)*SecDiffOperator;
        
        % Lineraization of (Dinternale_cxN05-thetaN05e*Dentropy_cxN05)*GxN05  part II
        Kmat3 = (-thetaAlgo*D2entropy_cxcxN05)*(GxN05(map.Voigt)*GxN05(map.Voigt)');
        
        % Material part of tangent operator
        ElasticityTensor = 2*Kmat2 + 4*Kmat3;
        TangentMaterialentr = 0.25*BN05'*ElasticityTensor*BN05;
        
        TangentMaterial = TangentMaterialentr + TangentMateriale;
        % Geometrical part of tangent operator
        GeometricalpartSave = dNx'*SAlgo*dNx;
        TangentGeometrical = zeros(numberOfMechanicalDOFs);
        for g = 1:DIM
            TangentGeometrical(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = GeometricalpartSave;
        end
        % Assembly of KXX
        kXX = 0.5*(TangentMaterial + TangentGeometrical);
        
        DQ_thetaN05 = -k0*(1/cxN05*GxN05*dNx*eye(size(dNx,2)));
        if flagDiscreteGradient(k,2) == false
            kXT = 0.5*kron(BN05'*(Dentropy_cxN05*De_theta*Dentropy_theta^(-2)*D2entropy_thetatheta)*GxN05(map.Voigt),N(k,:));
            kTT = N(k,:)'*1/DT*N(k,:) - 0.5*De_theta^(-1)*dNx'*DQ_thetaN05 - 0.5*N(k,:)'*Dentropy_theta^(-2)*D2entropy_thetatheta*(DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*Dentropy_cxN05)*N(k,:);
        else
            kXT = kron(BN05'*(Dentropy_cxN05*De_theta*Dentropy_theta^(-2)*D2entropy_thetatheta)*GxN05(map.Voigt),N(k,:));
            kTT = N(k,:)'*1/DT*N(k,:) - 0.5*De_theta^(-1)*dNx'*DQ_thetaN05-N(k,:)'*Dentropy_theta^(-2)*D2entropy_thetatheta*(DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*Dentropy_cxN05)*N(k,:);
        end
        
        % KThetaX
        kTX1 = N(k,:)'*(Dentropy_theta)^(-1)*DotCx(map.Voigt)'*(GxN05(map.Voigt)'*Isym)'*D2entropy_cxcxN05*(BN05'*GxN05(map.Voigt))';
        temp = wedge(CxN05,DotCx);
        kTX2 = N(k,:)'*(Dentropy_theta)^(-1)*Dentropy_cxN05*(BN05'*temp(map.Voigt))';
        kTX3 = N(k,:)'*(Dentropy_theta)^(-1)*Dentropy_cxN05*(1/DT*BN1'*GxN05(map.Voigt))';
        kTX4 = -dNx'*De_theta^(-1)*k0*cxN05^(-2)*(GxN05*(dNx*thetaN05'))*(BN05'*GxN05(map.Voigt))';
        Q_1G=zeros(3,3);
        Q_2G=zeros(3,3);
        Q_3G=zeros(3,3);
        I=eye(3);
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
        kTX5 = -(De_theta)^(-1)*dNx'*Q_Gv*Isym*0.5*SecDiffOperator*BN05;
        kTX = 0.5*kTX1 + 0.5*kTX2 + kTX3 + 0.5*kTX4 + 0.5*kTX5;
        %% summation of tangent
        KXX = KXX + kXX*detJ*wp(k);
        KXT = KXT + kXT*detJ*wp(k);
        KTX = KTX + kTX*detJ*wp(k);
        KTT = KTT + kTT*detJ*wp(k);
        
        ePot = [HelmholtzN1 EntropyN1 InternalEnergyN1 DeltaInternale];
        end
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