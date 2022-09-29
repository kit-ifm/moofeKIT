function [rData, kData, elementEnergy, array] = mixedSCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% MIXEDSCMOONEYRIVLINFULLCOUPLEDENDPOINT Element routine of class solidElectroClass.
%
% FORMULATION
% This is a 'mixed'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic  Mooney-Rivlin
% ('MooneyRivlin') model (nonlinear geometric and nonlinear material/
% stress-strain relation).
% The formulation is based on a polyconvexity inspired strain-energy
% density function.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
% In contrast to routine mixedSCMooneyRivlinEndpoint, this routine
% also considers thermoelectric coupling.
%
% CALL
% mixedSCMooneyRivlinEndpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [totalNumberOfFields,1] for residual data of
%        every field, here: (X, T, D, C, G, c, LambdaC, LambdaG, Lambdac)
% kData: cell-array of size [totalNumberOfFields, totalNumberOfFields] for
%        tangent data of every field, here: (X, T, D, C, G, c, LambdaC,
%        LambdaG, Lambdac)
% dofs: degrees of freedom (dofs) optionally manipulated data (numerical
%       tangent)
% array: structure for storage fe-data, for more information see
%        storageFEObject.initializeArrayStress
% stressTensor: structure for storage stress tensors (postprocessing), for
%               more information see storageFEObject.initializeArrayStress
% flagNumericalTangent: flag that indicates whether the function call
%                       happens during the computation of the numerical
%                       tangent or not.
%
% REFERENCE
% In preparation
%
% SEE ALSO
% mixedSCMooneyRivlinMidpoint
% mixedSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Felix Zaehringer, Marlon Franke

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
edAlphaN = obj.mixedFEObject.qN(e, :);

% internal DOFs
extractedDN = edAlphaN(1:3*numberOfInternalNodes).';
extractedCNv = edAlphaN(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
extractedGNv = edAlphaN(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
extractedcN = edAlphaN(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';
extractedCN1v = edAlphaN1(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
extractedGN1v = edAlphaN1(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
extractedcN1 = edAlphaN1(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
extractedLambdaCN1v = edAlphaN1(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
extractedLambdaGN1v = edAlphaN1(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
extractedLambdacN1 = edAlphaN1(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';

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
RCv = rData{4, 1};
RGv = rData{5, 1};
Rc = rData{6, 1};
RLambdaCv = rData{7, 1};
RLambdaGv = rData{8, 1};
RLambdac = rData{9, 1};
KXX=kData{1,1}; KXP=kData{1,2}; KXD=kData{1,3}; KXC=kData{1,4}; KXG=kData{1,5}; KXc=kData{1,6}; KXLambdaC=kData{1,7}; KXLambdaG=kData{1,8}; KXLambdac=kData{1,9}; 
KPX=kData{2,1}; KPP=kData{2,2}; KPD=kData{2,3}; KPC=kData{2,4}; KPG=kData{2,5}; KPc=kData{2,6}; KPLambdaC=kData{2,7}; KPLambdaG=kData{2,8}; KPLambdac=kData{2,9}; 
KDX=kData{3,1}; KDP=kData{3,2}; KDD=kData{3,3}; KDC=kData{3,4}; KDG=kData{3,5}; KDc=kData{3,6}; KDLambdaC=kData{3,7}; KDLambdaG=kData{3,8}; KDLambdac=kData{3,9}; 
KCX=kData{4,1}; KCP=kData{4,2}; KCD=kData{4,3}; KCC=kData{4,4}; KCG=kData{4,5}; KCc=kData{4,6}; KCLambdaC=kData{4,7}; KCLambdaG=kData{4,8}; KCLambdac=kData{4,9}; 
KGX=kData{5,1}; KGP=kData{5,2}; KGD=kData{5,3}; KGC=kData{5,4}; KGG=kData{5,5}; KGc=kData{5,6}; KGLambdaC=kData{5,7}; KGLambdaG=kData{5,8}; KGLambdac=kData{5,9}; 
KcX=kData{6,1}; KcP=kData{6,2}; KcD=kData{6,3}; KcC=kData{6,4}; KcG=kData{6,5}; Kcc=kData{6,6}; KcLambdaC=kData{6,7}; KcLambdaG=kData{6,8}; KcLambdac=kData{6,9}; 
KLambdaCX=kData{7,1}; KLambdaCP=kData{7,2}; KLambdaCD=kData{7,3}; KLambdaCC=kData{7,4}; KLambdaCG=kData{7,5}; KLambdaCc=kData{7,6}; KLambdaCLambdaC=kData{7,7}; KLambdaCLambdaG=kData{7,8}; KLambdaCLambdac=kData{7,9}; 
KLambdaGX=kData{8,1}; KLambdaGP=kData{8,2}; KLambdaGD=kData{8,3}; KLambdaGC=kData{8,4}; KLambdaGG=kData{8,5}; KLambdaGc=kData{8,6}; KLambdaGLambdaC=kData{8,7}; KLambdaGLambdaG=kData{8,8}; KLambdaGLambdac=kData{8,9}; 
KLambdacX=kData{9,1}; KLambdacP=kData{9,2}; KLambdacD=kData{9,3}; KLambdacC=kData{9,4}; KLambdacG=kData{9,5}; KLambdacc=kData{9,6}; KLambdacLambdaC=kData{9,7}; KLambdacLambdaG=kData{9,8}; KLambdacLambdac=kData{9,9};

%% DOFs
numberOfXDOFs = numel(RX);
numberOfDDOFs = numel(RD);

Isym = [eye(3) zeros(3); zeros(3) 2*eye(3)];

%% Gauss loop
for k = 1:numberOfGausspoints
    index = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);

    %% Independent strains/conjungated stresses
    DN = reshape(extractedDN, 3, []) * M(k, :)';
    CNv = reshape(extractedCNv, 6, []) * M(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    GNv = reshape(extractedGNv, 6, []) * M(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    cN = extractedcN.' * M(k, :)';
    DN1 = reshape(extractedDN1, 3, []) * M(k, :)';
    CN1v = reshape(extractedCN1v, 6, []) * M(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    cN1 = extractedcN1.' * M(k, :)';
    lambdaCN1v = reshape(extractedLambdaCN1v, 6, []) * M(k, :)';
    lambdaCN1 = voigtToMatrix(lambdaCN1v, 'stress');
    lambdaGN1v = reshape(extractedLambdaGN1v, 6, []) * M(k, :)';
    lambdaGN1 = voigtToMatrix(lambdaGN1v, 'stress');
    lambdacN1 = extractedLambdacN1.' * M(k, :)';

    %% Kinematic quantities
    FxN1 = edN1*dNx';
    CxN1 = FxN1'*FxN1;

    % Electrical quantities
    EN1 = -dNx*phiN1';

    % B-matrix mechanical (current configuration)
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

    % B-matrix electrical (current configuration)
    BEN1 = zeros(6,numberOfDDOFs);
    BEN1(1,1:3:end) = 2*DN1(1)*M(k,:);
    BEN1(2,2:3:end) = 2*DN1(2)*M(k,:);
    BEN1(3,3:3:end) = 2*DN1(3)*M(k,:);
    BEN1(4,1:3:end) = DN1(2)*M(k,:);
    BEN1(4,2:3:end) = DN1(1)*M(k,:);
    BEN1(5,2:3:end) = DN1(3)*M(k,:);
    BEN1(5,3:3:end) = DN1(2)*M(k,:);
    BEN1(6,1:3:end) = DN1(3)*M(k,:);
    BEN1(6,3:3:end) = DN1(1)*M(k,:);

    %% Strain energy function
    WmechN1 = a*(trace(CN1)-3) + b*(trace(GN1)-3);
    WvolN1 = c/2*(sqrt(cN1)-1)^2 - d*log(sqrt(cN1));
    WelectroMech = 1/(2*er*e0*cN1^(1/2))*DN1'*CxN1*DN1;
    % Energies
    elementEnergy.ePotMechanical  = elementEnergy.ePotMechanical + (WmechN1 + WvolN1)*detJ*gaussWeight(k);
    elementEnergy.ePotElectrical  = elementEnergy.ePotElectrical + (WelectroMech)*detJ*gaussWeight(k);
    elementEnergy.ePotD0timesGRADphi = elementEnergy.ePotD0timesGRADphi - DN1'*EN1*detJ*gaussWeight(k);

    %% Derivatives of the strain energy function
    DW_C = a*eye(3) + 1/(2*er*e0*cN1^(1/2))*(DN1*DN1');
    DW_G = b*eye(3);
    DW_c = c/2*(1-1/sqrt(cN1)) - d/(2*cN1) - 1/(4*er*e0)*cN1^(-3/2)*DN1'*(CN1*DN1);
    DW_D = CN1*DN1/(er*e0*cN1^(1/2));


    if ~computePostData

        %% Second derivatives of the strain energy function
        D2W_Cc = -1/(4*e0*er*cN1^(3/2))*(DN1*DN1');
        D2W_CD = 1/(2*er*e0*cN1^(1/2))*BEN1;
        
        D2W_cC = -1/(4*e0*er*cN1^(3/2))*DN1*DN1';
        D2W_cc = c/(4*cN1^(3/2))+d/(2*cN1^2)+3/(8*er*e0*cN1^(5/2))*DN1'*(CN1*DN1);
        D2W_cD = -1/(2*er*e0*cN1^(3/2))*CN1*DN1;
        
        D2W_DC = (1/(2*er*e0*cN1^(1/2))*BEN1)';
        D2W_Dc = -1/(2*er*e0*cN1^(3/2))*(CN1*DN1)';
        D2W_DD = 1/(er*e0*cN1^(1/2))*CN1;
        
        %% Residual densities
        % impuls evolution
        RX = RX + BN1'*lambdaCN1(map.Voigt)*detJ*gaussWeight(k);
        % electrical evolution
        % volume charge
        if isa(rho0,'function_handle')
            XGP = kron(N(k,:),eye(dimension))*edRef(:);
            rho0GP = rho0(XGP(1),XGP(2),XGP(3));
        else
            rho0GP = rho0;
        end
        RP = RP + (dNx'*DN1 + materialObject.timeFunctionRhoSource(time)*N(k,:)'*rho0GP)*detJ*gaussWeight(k);
        % constitutive compatibility conditions
        rC  = DW_C-lambdaCN1+wedge(lambdaGN1,CN1)+1/3*lambdacN1*GN1;
        RC = kron(M(k,:),Isym*rC(map.Voigt))*detJ*gaussWeight(k);
        RCv  = RCv + RC(:);
        rG  = DW_G-lambdaGN1+1/3*lambdacN1*CN1;
        RG = kron(M(k,:),Isym*rG(map.Voigt))*detJ*gaussWeight(k);
        RGv  = RGv + RG(:);
        Rc  = Rc + kron(M(k,:),(DW_c-lambdacN1))'*detJ*gaussWeight(k);
        % geometrically compatibility conditions
        rLambdaC = CxN1 - CN1;
        RLambdaC = kron(M(k,:),Isym*rLambdaC(map.Voigt))*detJ*gaussWeight(k);
        RLambdaCv = RLambdaCv + RLambdaC(:);
        rLambdaG = 0.5*wedge(CN1,CN1)-GN1;
        RLambdaG = kron(M(k,:),Isym*rLambdaG(map.Voigt))*detJ*gaussWeight(k);
        RLambdaGv = RLambdaGv + RLambdaG(:);
        RLambdac = RLambdac + kron(M(k,:),(1/3*sum(sum(GN1.*CN1))-cN1))'*detJ*gaussWeight(k);
        % electrical displacement evolution
        rD = kron(M(k,:),(DW_D-EN1))*detJ*gaussWeight(k);
        RD = RD + rD(:);
    
        %% Tangent
        kXX=zeros(numberOfXDOFs,numberOfXDOFs);
        A1 = 2*dNx'*lambdaCN1*dNx;
        for g = 1:dimension
            kXX(g:dimension:numberOfXDOFs,g:dimension:numberOfXDOFs) = A1;
        end
        KXX = KXX + kXX*detJ*gaussWeight(k);
        KXLambdaC = KXLambdaC + kron(M(k,:),BN1')*detJ*gaussWeight(k);
        KPD = KPD + kron(M(k,:),dNx')*detJ*gaussWeight(k);
        KDP = KDP + kron(M(k,:),dNx')'*detJ*gaussWeight(k);
        KDC = KDC + kron(M(k,:),(Isym*D2W_DC')')*detJ*gaussWeight(k);
        KDc = KDc + kron(M(k,:)'*M(k,:),D2W_Dc')*detJ*gaussWeight(k);
        KDD = KDD + kron(M(k,:)'*M(k,:),D2W_DD)*detJ*gaussWeight(k);
        KCLambdaC = KCLambdaC-kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        D=CN1;
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)	0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        KCLambdaG=KCLambdaG + kron(M(k,:)'*M(k,:),SecDiffOperator)*detJ*gaussWeight(k);
        D=lambdaGN1;
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)	0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        KCC= KCC + kron(M(k,:)'*M(k,:),SecDiffOperator)*detJ*gaussWeight(k);
        KCLambdac = KCLambdac + 1/3*kron(M(k,:)'*M(k,:),Isym*GN1(map.Voigt))*detJ*gaussWeight(k);
        KCG= KCG + 1/3*lambdacN1*kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        KCc = KCc + kron(M(k,:)'*M(k,:),D2W_Cc(map.Voigt))*detJ*gaussWeight(k);
        KCD = KCD + kron(M(k,:)',Isym*D2W_CD)*detJ*gaussWeight(k);
        KGLambdaG = KGLambdaG-kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        KGLambdac = KGLambdac + 1/3*kron(M(k,:)'*M(k,:),Isym*CN1(map.Voigt))*detJ*gaussWeight(k);
        KGC = KGC + 1/3*lambdacN1*kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        Kcc = Kcc + kron(M(k,:)'*M(k,:),1)*D2W_cc*detJ*gaussWeight(k);
        KcLambdac = KcLambdac -kron(M(k,:)'*M(k,:),1)*detJ*gaussWeight(k);
        KcC = KcC + kron(M(k,:)'*M(k,:),(D2W_cC(map.Voigt))')*detJ*gaussWeight(k);
        KcD = KcD + kron(M(k,:)'*M(k,:),D2W_cD)'*detJ*gaussWeight(k);
        KLambdaCX = KLambdaCX + kron(M(k,:)',BN1)*detJ*gaussWeight(k);
        KLambdaCC = KLambdaCC -kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        D=CN1;
        SecDiffOperator=[0          D(3,3)  D(2,2)      0        -D(3,2)  0      ;
            D(3,3)     0       D(1,1)      0        0        -D(3,1) ;
            D(2,2)     D(1,1)  0           -D(2,1)  0        0       ;
            0          0       -D(2,1)     -D(3,3)  D(3,1)   D(2,3)  ;
            -D(3,2)	0       0           D(3,1)   -D(1,1)  D(1,2)  ;
            0          -D(3,1)	0           D(3,2)   D(2,1)   -D(2,2) ];
        SecDiffOperator(1:3,4:6) = 2*SecDiffOperator(1:3,4:6);
        SecDiffOperator(4:6,1:6) = 2*SecDiffOperator(4:6,1:6);
        KLambdaGC=KLambdaGC + kron(M(k,:)'*M(k,:),SecDiffOperator)*detJ*gaussWeight(k);
        KLambdaGG=KLambdaGG -kron(M(k,:)'*M(k,:),Isym)*detJ*gaussWeight(k);
        KLambdacG=KLambdacG + 1/3*kron(M(k,:)'*M(k,:),(Isym*CN1(map.Voigt))')*detJ*gaussWeight(k);
        KLambdacC=KLambdacC + 1/3*kron(M(k,:)'*M(k,:),(Isym*GN1(map.Voigt))')*detJ*gaussWeight(k);
        KLambdacc =KLambdacc -kron(M(k,:)'*M(k,:),1)*detJ*gaussWeight(k);
    else

        %% Stress computation
        SN1 = 2 * lambdaCN1;
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
    rData{4} = RCv;
    rData{5} = RGv;
    rData{6} = Rc;
    rData{7} = RLambdaCv;
    rData{8} = RLambdaGv;
    rData{9} = RLambdac;
    kData{1, 1} = KXX;
    kData{1, 2} = KXP;
    kData{1, 3} = KXD;
    kData{1, 4} = KXC;
    kData{1, 5} = KXG;
    kData{1, 6} = KXc;
    kData{1, 7} = KXLambdaC;
    kData{1, 8} = KXLambdaG;
    kData{1, 9} = KXLambdac;
    kData{2, 1} = KPX;
    kData{2, 2} = KPP;
    kData{2, 3} = KPD;
    kData{2, 4} = KPC;
    kData{2, 5} = KPG;
    kData{2, 6} = KPc;
    kData{2, 7} = KPLambdaC;
    kData{2, 8} = KPLambdaG;
    kData{2, 9} = KPLambdac;
    kData{3, 1} = KDX;
    kData{3, 2} = KDP;
    kData{3, 3} = KDD;
    kData{3, 4} = KDC;
    kData{3, 5} = KDG;
    kData{3, 6} = KDc;
    kData{3, 7} = KDLambdaC;
    kData{3, 8} = KDLambdaG;
    kData{3, 9} = KDLambdac;
    kData{4, 1} = KCX;
    kData{4, 2} = KCP;
    kData{4, 3} = KCD;
    kData{4, 4} = KCC;
    kData{4, 5} = KCG;
    kData{4, 6} = KCc;
    kData{4, 7} = KCLambdaC;
    kData{4, 8} = KCLambdaG;
    kData{4, 9} = KCLambdac;
    kData{5, 1} = KGX;
    kData{5, 2} = KGP;
    kData{5, 3} = KGD;
    kData{5, 4} = KGC;
    kData{5, 5} = KGG;
    kData{5, 6} = KGc;
    kData{5, 7} = KGLambdaC;
    kData{5, 8} = KGLambdaG;
    kData{5, 9} = KGLambdac;
    kData{6, 1} = KcX;
    kData{6, 2} = KcP;
    kData{6, 3} = KcD;
    kData{6, 4} = KcC;
    kData{6, 5} = KcG;
    kData{6, 6} = Kcc;
    kData{6, 7} = KcLambdaC;
    kData{6, 8} = KcLambdaG;
    kData{6, 9} = KcLambdac;
    kData{7, 1} = KLambdaCX;
    kData{7, 2} = KLambdaCP;
    kData{7, 3} = KLambdaCD;
    kData{7, 4} = KLambdaCC;
    kData{7, 5} = KLambdaCG;
    kData{7, 6} = KLambdaCc;
    kData{7, 7} = KLambdaCLambdaC;
    kData{7, 8} = KLambdaCLambdaG;
    kData{7, 9} = KLambdaCLambdac;
    kData{8, 1} = KLambdaGX;
    kData{8, 2} = KLambdaGP;
    kData{8, 3} = KLambdaGD;
    kData{8, 4} = KLambdaGC;
    kData{8, 5} = KLambdaGG;
    kData{8, 6} = KLambdaGc;
    kData{8, 7} = KLambdaGLambdaC;
    kData{8, 8} = KLambdaGLambdaG;
    kData{8, 9} = KLambdaGLambdac;
    kData{9, 1} = KLambdacX;
    kData{9, 2} = KLambdacP;
    kData{9, 3} = KLambdacD;
    kData{9, 4} = KLambdacC;
    kData{9, 5} = KLambdacG;
    kData{9, 6} = KLambdacc;
    kData{9, 7} = KLambdacLambdaC;
    kData{9, 8} = KLambdacLambdaG;
    kData{9, 9} = KLambdacLambdac;
end
end