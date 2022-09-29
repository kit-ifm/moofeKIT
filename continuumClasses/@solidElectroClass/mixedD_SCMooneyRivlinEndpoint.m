function [rData, kData, elementEnergy, array] = mixedD_SCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof(e, :);
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
M = mixedFEObject.shapeFunctionObject.N;
% material data and voigt notation
numberOfInternalNodes = size(M, 2);
dimension = obj.dimension;
DT = setupObject.timeStepSize;
time = setupObject.time;

% map
map.Voigt = [1, 5, 9, 4, 8, 7]';
map.VoigtInv = [1, 4, 6; 4, 2, 5; 6, 5, 3];
map.VoigtFull = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
map.VoigtFullInv = [1, 4, 6; 7, 2, 5; 9, 8, 3];
map.Isym = [eye(3), zeros(3); zeros(3), 2 * eye(3)];
map.Isyminv = [eye(3), zeros(3); zeros(3), 1 / 2 * eye(3)];

% initialize energy
elementEnergy.ePotMechanical  = 0;
elementEnergy.ePotElectrical  = 0;
elementEnergy.ePotD0timesGRADphi = 0;

%% extract dofs for element
edN1 = dofs.edN1;
phiN1 = dofs.phiN1;
edAlphaN1 = dofs.edAlphaN1;

% internal DOFs
extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';

%
edRef = obj.qR(edof, 1:dimension)';

%% Jacobi-Matrix
J = edRef * dNr';
JN1 = edN1 * dNr';

%% Material parameters
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
e0 = materialObject.e0;
er = materialObject.e1;
rho0 = materialObject.rhoSource;

%% Initialization
RX = rData{1, 1};
RP = rData{2, 1};
RD = rData{3, 1};
KXX=kData{1,1}; KXP=kData{1,2}; KXD=kData{1,3};
KPX=kData{2,1}; KPP=kData{2,2}; KPD=kData{2,3};
KDX=kData{3,1}; KDP=kData{3,2}; KDD=kData{3,3};

%% DOFs
numberOfXDOFs = numel(RX);
numberOfPDOFs = numel(RP);

%% Gauss loop
for k = 1:numberOfGausspoints
    index = dimension*k-(dimension-1):dimension*k;
    detJ = det(J(:,index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,index)')\dNr(index,:);

    % Independent electrical displacement field
%     DN1 = zeros(3,1);
%     for m = 1:numberOfInternalNodes
%         DN1 = DN1 + M(k,m)*DN1e(e,ones(1,3)*(m - 1)*3 + (1:3))';
%     end
    DN1 = reshape(extractedDN1, 3, []) * M(k, :)';

    % Kinematic quantities
    FxN1 = edN1*dNx';
    CxN1 = FxN1'*FxN1;
    GxN1 = 0.5*wedge(CxN1,CxN1);
    cxN1 = det(CxN1);

    % Electrical quantities
    EN1 = -dNx*phiN1';

    % B-matrix (current configuration)
    BN1 = zeros(6,numberOfXDOFs);
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
    WmechN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3);
    WvolN1 = c/2*(sqrt(cxN1)-1)^2 - d*log(sqrt(cxN1));
    WelectroMech = 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;
    WN1 = WmechN1 + WvolN1 + WelectroMech;

    % energies
    elementEnergy.ePotMechanical  = elementEnergy.ePotMechanical + (WmechN1 + WvolN1)*detJ*gaussWeight(k);
    elementEnergy.ePotElectrical  = elementEnergy.ePotElectrical + (WelectroMech)*detJ*gaussWeight(k);
    elementEnergy.ePotD0timesGRADphi = elementEnergy.ePotD0timesGRADphi - DN1'*EN1*detJ*gaussWeight(k);


    %% Residual
    % Impuls evolution
    % Derivative of the strain energy function
    DW_CN1 = a*eye(3) + 1/(2*er*e0*cxN1^(1/2))*DN1*DN1';
    DW_GN1 = b*eye(3);
    DW_cN1 = c/2*(1-1/sqrt(cxN1)) - d/(2*cxN1) - 1/(4*er*e0)*cxN1^(-3/2)*DN1'*(CxN1*DN1);
    DW_DN1 = CxN1*DN1/(er*e0*cxN1^(1/2));
    % Conjungated stresses
    SN1 = 2*(DW_CN1 + wedge(DW_GN1,CxN1) + DW_cN1*GxN1);

    if ~computePostData

        %% Residual densities
        % impuls evolution
        RX = RX + BN1'*0.5*SN1(map.Voigt)*detJ*gaussWeight(k);
        % volume charge
        if isa(rho0,'function_handle')
            XGP = kron(N(k,:),eye(dimension))*edRef(:);
            rho0GP = rho0(XGP(1),XGP(2),XGP(3));
        else
            rho0GP = rho0;
        end
        % electrical evolution
        RP = RP + (dNx'*DN1 + materialObject.timeFunctionRhoSource(time)*N(k,:)'*rho0GP)*detJ*gaussWeight(k);
        % electrical displacement evolution
        rD = kron(M(k,:),(DW_DN1-EN1));
        RD = RD + rD(:)*detJ*gaussWeight(k);

        %% Tangent
        % Lineraization of DW_CN1
        D2W_CN1cN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        Kmat0 = D2W_CN1cN1(map.Voigt)*GxN1(map.Voigt)';
        % Lineraization of wedge(DW_GN1,CxN1)
        D = DW_GN1;
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0  ;
            D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
            D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
            0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
            -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
            0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
        SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
        SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
        Kmat1 = SecDiffOperator;
        % Lineraization of DW_cN1*GxN1  part I
        D = (CxN1);
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0   ;
            D(3,3)     0       D(1,1)      0         0       -D(3,1)  ;
            D(2,2)     D(1,1)  0           -D(2,1)	  0       0        ;
            0          0       -D(2,1)     -D(3,3)   D(3,1)  D(2,3)   ;
            -D(3,2)	    0       0           D(3,1)   -D(1,1)  D(1,2)   ;
            0          -D(3,1)	0           D(3,2)    D(2,1)  -D(2,2)] ;
        SecDiffOperator(1:3,1:6) = 2*SecDiffOperator(1:3,1:6);
        SecDiffOperator(4:6,1:3) = 2*SecDiffOperator(4:6,1:3);
        Kmat2 = DW_cN1*SecDiffOperator;
        % Lineraization of DW_cN1*GxN1  part II
        % Second derivative of the strain energy function
        D2W_cc = c/(4*cxN1^(3/2)) + d/(2*cxN1^2) + 3/(8*er*e0)*cxN1^(-5/2)*DN1'*(CxN1*DN1);
        Kmat3 = D2W_cc*(GxN1(map.Voigt)*GxN1(map.Voigt)');
        D2W_cN1CN1 = -1/(4*er*e0)*cxN1^(-3/2)*DN1*DN1';
        Kmat4 = GxN1(map.Voigt)*D2W_cN1CN1(map.Voigt)';
        % Material part of tangent operator
        ElasticityTensor = 4*Kmat0 + 2*Kmat1  + 2*Kmat2 + 4*Kmat3 + 4*Kmat4;
        TangentMaterial = 0.25*BN1'*ElasticityTensor*BN1;
        % Geometrical part of tangent operator
        GeometricalpartSave = dNx'*SN1*dNx;
        TangentGeometrical = zeros(numberOfXDOFs);
        for g = 1:dimension
            TangentGeometrical(g:dimension:numberOfXDOFs,g:dimension:numberOfXDOFs) = GeometricalpartSave;
        end
        KXX = KXX + (TangentMaterial + TangentGeometrical)*detJ*gaussWeight(k);
        % KXP
        KXP = zeros(numberOfXDOFs,numberOfPDOFs);
        % KXD
        BEN1 = zeros(6,3*numberOfInternalNodes);
        BEN1(1,1:3:end) = 2*DN1(1)*M(k,:);
        BEN1(2,2:3:end) = 2*DN1(2)*M(k,:);
        BEN1(3,3:3:end) = 2*DN1(3)*M(k,:);
        BEN1(4,1:3:end) = DN1(2)*M(k,:);
        BEN1(4,2:3:end) = DN1(1)*M(k,:);
        BEN1(5,2:3:end) = DN1(3)*M(k,:);
        BEN1(5,3:3:end) = DN1(2)*M(k,:);
        BEN1(6,1:3:end) = DN1(3)*M(k,:);
        BEN1(6,3:3:end) = DN1(1)*M(k,:);
        D2W_cN1DN1 = -1/(2*er*e0)*cxN1^(-3/2)*(CxN1*DN1);
        D2W_DN1cN1 = -CxN1*DN1/(2*er*e0)*cxN1^(-3/2);
        temp1 = BN1'*GxN1(map.Voigt);
        temp2 = kron(M(k,:),D2W_DN1cN1);
        KXD = KXD + (BN1'*0.5*2*1/(2*er*e0*cxN1^(1/2))*BEN1 + temp1*temp2(:)')*detJ*gaussWeight(k);
        % KPX
        KPX = zeros(numberOfPDOFs,numberOfXDOFs);
        % KPP
        KPP = zeros(numberOfPDOFs,numberOfPDOFs);
        % KPD
        KPD = KPD + kron(M(k,:),dNx')*detJ*gaussWeight(k);
        % KDX
        KDX = KDX + (BEN1'*0.5*2*1/(2*er*e0*cxN1^(1/2))*BN1 + temp2(:)*temp1')*detJ*gaussWeight(k);
        % KDP
        KDP = KDP + kron(M(k,:),dNx')'*detJ*gaussWeight(k);
        % KDD
        KDD = KDD + kron(M(k,:)'*M(k,:),CxN1/(er*e0*cxN1^(1/2)))*detJ*gaussWeight(k);
    else

        %% Stress computation
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        stressTensor.D = DN1;
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    rData{2} = RP;
    rData{3} = RD;
    kData{1, 1} = KXX;
    kData{1, 2} = KXP;
    kData{1, 3} = KXD;
    kData{2, 1} = KPX;
    kData{2, 2} = KPP;
    kData{2, 3} = KPD;
    kData{3, 1} = KDX;
    kData{3, 2} = KDP;
    kData{3, 3} = KDD;
end
end
