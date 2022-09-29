function [rData, kData, elementEnergy, array] = mixedD_SCMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% 05.04.2019 M.F.

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;

% aquire general data
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
if strcmpi(mixedFEObject.typeShapeFunction,'sameOrder')
    N_D = mixedFEObject.shapeFunctionObject.N;
elseif strcmpi(mixedFEObject.typeShapeFunction,'detailedOrder')
    N_D = mixedFEObject.shapeFunctionObject.N{1};
end


numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);
numberOfXDOFs = numel(rData{1});
numberOfInternalNodes = size(N_D, 2);

dimension = obj.dimension;

DT = setupObject.timeStepSize;

I6 = [eye(3), zeros(3); zeros(3), 2 * eye(3)];

% aquire material data
a = materialObject.a;
b = materialObject.b;
c1 = materialObject.c1;
d1 = materialObject.d1;
c2 = materialObject.c2;
d2 = materialObject.d2;
e0 = materialObject.e0;
er = materialObject.e1;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;
rho0 = materialObject.rhoSource;
R = materialObject.RSource;

% aquire the nodal values of the variables for the current element
edRef = obj.qR(edof, 1:dimension).';

edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 1 / 2 * (edN + edN1);

phiN = obj.qN(edof, dimension+1).';
phiN1 = dofs.phiN1;
phiN05 = 0.5 * (phiN + phiN1);

thetaN = obj.qN(edof, dimension+2).';
thetaN1 = dofs.thetaN1;
thetaN05 = 0.5 * (thetaN1 + thetaN);

edAlphaN = obj.mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedDN = edAlphaN(1:3*numberOfInternalNodes);
extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes);

% Jacobian matrices
J = edRef * dNr';
JN1 = edN1 * dNr';

% initialize residual
RX = rData{1, 1};
RP = rData{2, 1};
RT = rData{3, 1};
RD = rData{4, 1};

% initialize tangent
KXX = kData{1, 1};
KXT = kData{1, 3};
KXD = kData{1, 4};
KPD = kData{2, 4};
KTX = kData{3, 1};
KTT = kData{3, 3};
KDX = kData{4, 1};
KDP = kData{4, 2};
KDD = kData{4, 4};

% initialize elementEnergy
elementEnergy.internalEnergy = 0;
elementEnergy.helmholtz = 0;
elementEnergy.DE = 0;
elementEnergy.TS = 0;
elementEnergy.S = 0;
elementEnergy.deltaS = 0;
elementEnergy.deltaU = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 4)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    index = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, index)');
    detJ1 = det(JN1(:, index)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);
    % Temperature
    thetaN1e = N(k, :) * thetaN1.';
    thetaNe = N(k, :) * thetaN.';
    thetaN05e = 0.5 * (thetaN1e+thetaNe);
    DotTheta = (thetaN1e-thetaNe) / DT;
    % Deformation Gradient
    FxN1 = edN1 * dNx';
    FxN = edN * dNx';
    FxN05 = edN05 * dNx';
    CxN1 = FxN1.' * FxN1;
    CxN = FxN.' * FxN;
    CxN05 = FxN05.' * FxN05;
    CxN05Algo = 0.5 * (CxN + CxN1);
    DotCx = (CxN1 - CxN) / DT;
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    GxN = 0.5 * wedge(CxN, CxN);
    GxN05 = 0.5 * wedge(CxN05, CxN05);
    GxN05Algo = 0.5 * (GxN + GxN1);
    cxN1 = det(CxN1);
    cxN = det(CxN);
    cxN05 = det(CxN05);
    % Electrical quantities
    EN1 = -dNx * phiN1.';
    EN = -dNx * phiN.';
    EN05 = -dNx * phiN05.';
    DN1 = reshape(extractedDN1, 3, []) * N_D(k, :)';
    DN = reshape(extractedDN, 3, []) * N_D(k, :)';
    DN05 = 0.5 * (DN + DN1);

    % B-matrix
    BN1 = BMatrix(dNx, FxN1);
    BN05 = BMatrix(dNx, FxN05);

    % ENERGY
    PsiTN = kappa * (thetaNe - thetaR - thetaNe * log(thetaNe/thetaR));
    PsiTMN = -dimension * beta * c2 * (thetaNe - thetaR) * (cxN - 1);
    PsiEMN = 1 / (2 * er * e0 * cxN^(1 / 2)) * DN' * CxN * DN;
    PsiMIsoN = a * (trace(CxN) - 3) + b * (trace(GxN) - 3);
    PsiMVolN = c1 / 2 * (sqrt(cxN) - 1)^2 - d1 * log(sqrt(cxN));
    PsiN = PsiTN + PsiTMN + PsiEMN + PsiMIsoN + PsiMVolN;
    s0N = mooneyRivlin(obj, CxN, GxN, cxN, DN, thetaNe, 3);
    u0N = PsiN - DN' * EN + thetaNe * s0N;
    PsiTN1 = kappa * (thetaN1e-thetaR - thetaN1e * log(thetaN1e/thetaR));
    %     PsiTMN1 = DIM*beta*c2*(thetaN1e - thetaR)*(cxN1 - 1);
    PsiTMN1 = -dimension * beta * c2 * (thetaN1e-thetaR) * (cxN1 - 1);
    PsiEMN1 = 1 / (2 * er * e0 * cxN1^(1 / 2)) * DN1' * CxN1 * DN1;
    PsiMIsoN1 = a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3);
    PsiMVolN1 = c1 / 2 * (sqrt(cxN1) - 1)^2 - d1 * log(sqrt(cxN1));
    PsiN1 = PsiTN1 + PsiTMN1 + PsiEMN1 + PsiMIsoN1 + PsiMVolN1;
    s0N1 = mooneyRivlin(obj, CxN1, GxN1, cxN1, DN1, thetaN1e, 3);
    u0N1 = PsiN1 - DN1' * EN1 + thetaN1e * s0N1;
    %     MR-Mike
    %     Psi = a*(trace(CN1)-DIM) + b*(1/2*((trace(CN1))^2 - trace((CN1)^2))-DIM) - d1*log(sqrt(cxN1)) + kappa*(thetaN1 - thetaR -thetaN1*log(thetaN1/thetaR)) + DIM*beta*(thetaN1 - thetaR)*d1*I3N1^(-1);
    %     MR-ThermoMech
    %     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR)) - DIM*beta*(thetaN1e - thetaR)*(c2*(sqrt(cxN1) - 1) - d2/sqrt(cxN1));
    %     MR-ThermoElectroMech
    %     PsiN1 = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) + c1/2*(sqrt(cxN1)-1)^2 - d1*log(sqrt(cxN1)) + kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR)) - DIM*beta*(thetaN1e - thetaR)*c2*(cxN1 - 1) + 1/(2*er*e0*cxN1^(1/2))*DN1'*CxN1*DN1;

    elementEnergy.internalEnergy = elementEnergy.internalEnergy + u0N1 * detJ * gaussWeight(k);
    elementEnergy.helmholtz = elementEnergy.helmholtz + PsiN1 * detJ * gaussWeight(k);
    elementEnergy.DE = elementEnergy.DE + DN1' * EN1 * detJ * gaussWeight(k);
    elementEnergy.TS = elementEnergy.TS + thetaN1e * s0N1 * detJ * gaussWeight(k);
    elementEnergy.S = elementEnergy.S + s0N1 * detJ * gaussWeight(k);
    elementEnergy.deltaS = elementEnergy.deltaS + (s0N1 - s0N) * detJ * gaussWeight(k);
    elementEnergy.deltaU = elementEnergy.deltaU + (u0N1 - u0N) * detJ * gaussWeight(k);

    % DERIVATIVES
    deltaC = CxN1 - CxN;
    normDeltaC = sqrt(innerProduct(deltaC, deltaC));
    deltac = cxN1 - cxN;
    normDeltac = abs(deltac);
    deltaD = DN1 - DN;
    normDeltaD = sqrt(deltaD'*deltaD);
    deltaT = thetaN1e-thetaNe;
    normDeltaT = abs(deltaT);

    if ~flagNumericalTangent
        discreteGradientCondition_C = (normDeltaC > 10 * eps);
        discreteGradientCondition_c = (normDeltac > 10 * eps);
        discreteGradientCondition_D = (normDeltaD > 10 * eps);
        discreteGradientCondition_T = (normDeltaT > 10 * eps);
    else
        discreteGradientCondition_C = flagDiscreteGradient(k, 1);
        discreteGradientCondition_c = flagDiscreteGradient(k, 2);
        discreteGradientCondition_D = flagDiscreteGradient(k, 3);
        discreteGradientCondition_T = flagDiscreteGradient(k, 4);
    end

    if discreteGradientCondition_C
        flagDiscreteGradient(k, 1) = 1;
        DW_C = a * eye(3) + 1 / (4 * er * e0) * (1 / sqrt(cxN1) * (DN1 * DN1.') + 1 / sqrt(cxN) * (DN * DN.'));
    else
        flagDiscreteGradient(k, 1) = 0;
        DW_C = a * eye(3) + 1 / (2 * er * e0 * sqrt(cxN05)) * (DN05 * DN05.');
    end

    DW_G = b * eye(3);

    if discreteGradientCondition_c
        flagDiscreteGradient(k, 2) = 1;
        DW_c = d1 / deltac * log(sqrt(cxN/cxN1)) ...
            +c1 / (2 * deltac) * ((sqrt(cxN1) - 1)^2 - (sqrt(cxN) - 1)^2) ...
            +1 / (4 * er * e0 * deltac) * (1 / sqrt(cxN1) - 1 / sqrt(cxN)) * (DN.' * (CxN1 * DN) + DN1.' * (CxN * DN1)) ...
            -3 / 2 * beta * c2 * (thetaN1e+thetaNe - 2 * thetaR);
    else
        flagDiscreteGradient(k, 2) = 0;
        DW_c = -d1 / (2 * cxN05) + c1 / 2 * (1 - cxN05^(-1 / 2)) - 1 / (4 * er * e0 * cxN05^(3 / 2)) * DN05.' * CxN05 * DN05 - dimension * beta * c2 * (thetaN05e-thetaR);
    end

    if discreteGradientCondition_D
        flagDiscreteGradient(k, 3) = 1;
        DW_D = 1 / (2 * er * e0) * (CxN1 / sqrt(cxN1) + CxN / sqrt(cxN)) * DN05;
    else
        flagDiscreteGradient(k, 3) = 0;
        DW_D = 1 / (er * e0 * sqrt(cxN05)) * CxN05 * DN05;
    end

    if discreteGradientCondition_T
        flagDiscreteGradient(k, 4) = 1;
        DW_theta = -3 / 2 * beta * c2 * (cxN1 + cxN - 2) ...
            +kappa / deltaT * (thetaN1e * (1 - log(thetaN1e/thetaR)) - thetaNe * (1 - log(thetaNe/thetaR)));
    else
        flagDiscreteGradient(k, 4) = 0;
        DW_theta = -dimension * beta * c2 * (cxN05 - 1) - kappa * log(thetaN05e/thetaR);
    end

    % RESIDUAL
    CAlgo = CxN05Algo;
    GAlgo = 1 / 3 * (wedge(CxN05Algo, CxN05Algo) + GxN05Algo);
    etaAlgo = -DW_theta;
    EAlgo = DW_D;
    SAlgo = 2 * (DW_C + wedge(DW_G, CAlgo) + DW_c * GAlgo);
    QxN05 = -k0 * cxN05^(-1) * GxN05 * (dNx * thetaN05.');

    RX = RX + BN05.' * matrixToVoigt(SAlgo, 'stress') * detJ * gaussWeight(k);
    RP = RP + (dNx' * DN05) * detJ * gaussWeight(k);
    RD = RD + kron(N_D(k, :)', (EAlgo - EN05)) * detJ * gaussWeight(k);
    RT = RT + (N(k, :)' * (thetaN1e * s0N1 - thetaNe * s0N) / DT - N(k, :)' * (thetaN1e-thetaNe) / DT * etaAlgo - dNx' * QxN05) * detJ * gaussWeight(k);

    % TANGENT
    % KXX
    geometricalTangentPart = 1 / 2 * dNx' * SAlgo * dNx;
    geometricalTangent = zeros(numberOfXDOFs);
    for g = 1:dimension
        geometricalTangent(g:dimension:numberOfXDOFs, g:dimension:numberOfXDOFs) = geometricalTangentPart;
    end

    if flagDiscreteGradient(k, 1)
        D2W_C_cxN1 = -1 / (8 * er * e0 * cxN1^(3 / 2)) * (DN1 * DN1');
        kmat1 = BN05' * 2 * matrixToVoigt(D2W_C_cxN1, 'stress') * matrixToVoigt(GxN1, 'stress')' * 2 * BN1;
    else
        D2W_C_cxN05 = -1 / (4 * er * e0 * cxN05^(3 / 2)) * (DN05 * DN05');
        kmat1 = BN05' * 2 * matrixToVoigt(D2W_C_cxN05, 'stress') * matrixToVoigt(GxN05, 'stress')' * BN05;
    end
    kmat2 = BN05' * 2 * secDiffOperator(DW_G) * BN1;
    if flagDiscreteGradient(k, 2)
        D2W_c_cxN1 = -1 / 2 * d1 / (cxN1 * deltac) ...
            -d1 / (deltac)^2 * log(sqrt(cxN/cxN1)) ...
            +1 / 2 * c1 / deltac * (1 - 1 / sqrt(cxN1)) ...
            -1 / 2 * c1 / (deltac)^2 * ((sqrt(cxN1) - 1)^2 - (sqrt(cxN) - 1)^2) ...
            -1 / 8 * 1 / (er * e0 * (cxN1)^(3 / 2) * deltac) * (DN.' * (CxN1 * DN) + DN1.' * (CxN * DN1)) ...
            -1 / 4 * 1 / (er * e0 * deltac^2) * (1 / sqrt(cxN1) - 1 / sqrt(cxN)) * (DN.' * (CxN1 * DN) + DN1.' * (CxN * DN1));

        kmat3 = BN05' * 2 * matrixToVoigt(GAlgo, 'stress') * 1 / (4 * er * e0 * deltac) * (1 / sqrt(cxN1) - 1 / sqrt(cxN)) * DN' * BErs(DN) * 2 * BN1 ...
            +BN05' * 2 * D2W_c_cxN1 * (matrixToVoigt(GAlgo, 'stress') * matrixToVoigt(GxN1, 'stress')') * 2 * BN1;
    else
        D2W_c_cxN05 = d1 / (2 * cxN05^2) + c1 / (4 * cxN05^(3 / 2)) + 3 / (8 * er * e0 * cxN05^(5 / 2)) * DN05.' * CxN05 * DN05;

        kmat3 = -BN05' * 2 * matrixToVoigt(GAlgo, 'stress') * 1 / (4 * er * e0 * cxN05^(3 / 2)) * DN05' * BErs(DN05) * BN05 ...
            +BN05' * 2 * D2W_c_cxN05 * (matrixToVoigt(GAlgo, 'stress') * matrixToVoigt(GxN05, 'stress')') * BN05;
    end
    kmat4 = BN05' * 2 * DW_c * secDiffOperator(1/3*(CAlgo + 1 / 2 * CxN1)) * 2 * BN1;

    KXX = KXX + (geometricalTangent + kmat1 + kmat2 + kmat3 + kmat4) * detJ * gaussWeight(k);

    % KXD
    if flagDiscreteGradient(k, 2)
        D2W_c_DN1 = 1 / (2 * er * e0 * deltac) * (1 / sqrt(cxN1) - 1 / sqrt(cxN)) * CxN * DN1;
        kXD1 = kron(N_D(k, :), BN05'*1/(2 * er * e0 * sqrt(cxN1))*AErs(DN1));
        kXD2 = kron(N_D(k, :), BN05'*2*matrixToVoigt(GAlgo, 'stress')*D2W_c_DN1');
    else
        D2W_c_DN05 = -1 / (2 * er * e0 * cxN05^(3 / 2)) * CxN05 * DN05;
        kXD1 = kron(N_D(k, :), BN05'*2*1/(2 * er * e0 * sqrt(cxN05))*AErs(DN05)*0.5);
        kXD2 = kron(N_D(k, :), BN05'*2*matrixToVoigt(GAlgo, 'stress')*D2W_c_DN05'*0.5);
    end
    KXD = KXD + (kXD1 + kXD2) * detJ * gaussWeight(k);

    % KXT
    if flagDiscreteGradient(k, 2)
        D2W_c_thetaN1 = -dimension / 2 * beta * c2;
        kXT1 = kron(N(k, :), BN05'*2*D2W_c_thetaN1*matrixToVoigt(GAlgo, 'stress'));
    else
        D2W_c_thetaN05 = -dimension * beta * c2;
        kXT1 = kron(N(k, :), BN05'*2*D2W_c_thetaN05*matrixToVoigt(GAlgo, 'stress')*0.5);
    end
    KXT = KXT + kXT1 * detJ * gaussWeight(k);

    % =====

    % KPD
    KPD = KPD + 1 / 2 * kron(N_D(k, :), dNx') * detJ * gaussWeight(k);

    % =====

    % KDX
    if flagDiscreteGradient(k, 3)
        D2W_D_cN1 = -1 / (4 * er * e0 * cxN1^(3 / 2)) * CxN1 * DN05;
        kDX1 = kron(N_D(k, :)', D2W_D_cN1*matrixToVoigt(GxN1, 'stress')'*2*BN1);
        kDX2 = kron(N_D(k, :)', 1/(2 * er * e0 * sqrt(cxN1))*BErs(DN05)*2*BN1);
    else
        D2W_D_cN05 = -1 / (2 * er * e0 * (cxN05)^(3 / 2)) * CxN05 * DN05;
        kDX1 = kron(N_D(k, :)', D2W_D_cN05*matrixToVoigt(GxN05, 'stress')'*BN05);
        kDX2 = kron(N_D(k, :)', 1/(er * e0 * sqrt(cxN05))*BErs(DN05)*BN05);
    end
    KDX = KDX + (kDX1 + kDX2) * detJ * gaussWeight(k);

    % KDP
    KDP = KDP + 1 / 2 * kron(N_D(k, :)', dNx) * detJ * gaussWeight(k);

    % KDD
    if flagDiscreteGradient(k, 3)
        D2W_D_DN1 = 1 / (4 * er * e0) * (CxN1 / sqrt(cxN1) + CxN / sqrt(cxN));
        kDD1 = kron(N_D(k, :)'*N_D(k, :), D2W_D_DN1);
    else
        D2W_D_DN05 = 1 / (er * e0 * sqrt(cxN05)) * CxN05;
        kDD1 = kron(N_D(k, :)'*N_D(k, :), D2W_D_DN05*0.5);
    end
    KDD = KDD + kDD1 * detJ * gaussWeight(k);

    % =====

    % KTX
    DEtaN1_cxN1 = dimension * beta * c2;
    kTX1 = N(k, :)' * thetaN1e / DT * DEtaN1_cxN1 * matrixToVoigt(GxN1, 'stress')' * 2 * BN1;
    if flagDiscreteGradient(k, 4)
        D2W_theta_cxN1 = -dimension / 2 * beta * c2;
        DEtaAlgo_cxN1 = -D2W_theta_cxN1;
        kTX2 = -N(k, :)' * (thetaN1e-thetaNe) / DT * DEtaAlgo_cxN1 * matrixToVoigt(GxN1, 'stress')' * 2 * BN1;
    else
        D2W_theta_cxN05 = -dimension * beta * c2;
        DEtaAlgo_cxN05 = -D2W_theta_cxN05;
        kTX2 = -N(k, :)' * (thetaN1e-thetaNe) / DT * DEtaAlgo_cxN05 * matrixToVoigt(GxN05, 'stress')' * BN05;
    end
    kTX3 = -dNx' * k0 * cxN05^(-2) * GxN05 * (dNx * thetaN05') * matrixToVoigt(GxN05, 'stress')' * BN05;
    kTX4 = dNx' * k0 * cxN05^(-1) * BErs(dNx*thetaN05') * secDiffOperator2(CxN05) * BN05;
    KTX = KTX + (kTX1 + kTX2 + kTX3 + kTX4) * detJ * gaussWeight(k);

    % KTT
    DEtaN1_thetaN1 = kappa / thetaN1e;
    kTT1 = N(k, :)' * (s0N1 + thetaN1e * DEtaN1_thetaN1) / DT * N(k, :);
    if flagDiscreteGradient(k, 4)
        D2W_theta_thetaN1e = -kappa / deltaT * log(thetaN1e/thetaR) ...
            -kappa / deltaT^2 * (thetaN1e * (1 - log(thetaN1e/thetaR)) - thetaNe * (1 - log(thetaNe/thetaR)));
        DEtaAlgo_thetaN1 = -D2W_theta_thetaN1e;
        kTT2 = -N(k, :)' * (etaAlgo + (thetaN1e-thetaNe) * DEtaAlgo_thetaN1) / DT * N(k, :);
    else
        D2W_theta_thetaN05e = -kappa / thetaN05e;
        DEtaAlgo_thetaN05 = -D2W_theta_thetaN05e;
        kTT2 = -N(k, :)' * (etaAlgo + (thetaN1e-thetaNe) * DEtaAlgo_thetaN05 * 0.5) / DT * N(k, :);
    end
    kTT3 = dNx' * k0 * cxN05^(-1) * GxN05 * 0.5 * dNx;
    KTT = KTT + (kTT1 + kTT2 + kTT3) * detJ * gaussWeight(k);
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RP;
    rData{3} = RT;
    rData{4} = RD;

    % pass tangent
    kData{1, 1} = KXX;
    kData{1, 3} = KXT;
    kData{1, 4} = KXD;
    kData{2, 4} = KPD;
    kData{3, 1} = KTX;
    kData{3, 3} = KTT;
    kData{4, 1} = KDX;
    kData{4, 2} = KDP;
    kData{4, 4} = KDD;

    % pass flagDiscreteGradient
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
end
end


function out = mooneyRivlin(obj, C, G, c, D, theta, choice)
a = obj.materialObject.a;
b = obj.materialObject.b;
c1 = obj.materialObject.c1;
d1 = obj.materialObject.d1;
c2 = obj.materialObject.c2;
d2 = obj.materialObject.d2;
e0 = obj.materialObject.e0;
er = obj.materialObject.e1;
kappa = obj.materialObject.kappa;
beta = obj.materialObject.beta;
thetaR = obj.materialObject.thetaR;
dimension = obj.dimension;
%
PsiT = kappa * (theta - thetaR - theta * log(theta/thetaR));
PsiTM = -dimension * beta * c2 * (theta - thetaR) * (c - 1);
PsiEM = 1 / (2 * er * e0 * c^(1 / 2)) * D' * C * D;
PsiMIso = a * (trace(C) - 3) + b * (trace(G) - 3);
PsiMVol = c1 / 2 * (sqrt(c) - 1)^2 - d1 * log(sqrt(c));
Psi = PsiT + PsiTM + PsiEM + PsiMIso + PsiMVol;
%
DPsi_C = 1 / (2 * er * e0) * c^(-1 / 2) * D * D' + a * eye(3);
DPsi_c = -dimension * beta * c2 * (theta - thetaR) - 1 / (4 * er * e0 * c^(3 / 2)) * D' * C * D + c1 / 2 * (1 - c^(-1 / 2)) - d1 / (2 * c);
DPsi_D = 1 / (er * e0 * c^(1 / 2)) * C * D;
DPsi_t = -kappa * log(theta/thetaR) - dimension * beta * c2 * (c - 1);
%
D2Psi_Cc = -1 / (4 * er * e0) * c^(-3 / 2) * D * D';
%
s0 = -DPsi_t;
Ds0_c = dimension * beta * c2;
%
% u0 = Psi + theta*s0;
% u0 = Psi + theta*s0 - D'*E;
u0 = kappa * (theta - thetaR) + dimension * beta * c2 * thetaR * (c - 1) + 1 / (2 * er * e0 * c^(1 / 2)) * D' * C * D + a * (trace(C) - 3) + b * trace(G) - d1 * log(c^(1 / 2)) + c1 / 2 * (c^(1 / 2) - 1)^2;
Du0_C = 1 / (2 * er * e0 * c^(1 / 2)) * D * D' + a * eye(3);
Du0_c = dimension * beta * c2 * thetaR - 1 / (4 * er * e0 * c^(3 / 2)) * D' * C * D + c1 / 2 * (1 - c^(-1 / 2)) - d1 / (2 * c);
Du0_D = 1 / (er * e0 * c^(1 / 2)) * C * D;
%
if choice == 1
    out = Psi;
elseif choice == 2
    out = u0;
elseif choice == 3
    out = s0;
elseif choice == 4
    out = Du0_C;
elseif choice == 5
    out = Du0_c;
elseif choice == 6
    out = Du0_D;
elseif choice == 14
    out = DPsi_C;
elseif choice == 15
    out = DPsi_c;
elseif choice == 16
    out = DPsi_D;
elseif choice == 17
    out = DPsi_t;
elseif choice == 21
    out = D2Psi_Cc;
end
end

function out = AErs(a)
out = zeros(6, 3);
if isvector(a)
    out = [2 * a(1), 0, 0; 0, 2 * a(2), 0; 0, 0, 2 * a(3); a(2), a(1), 0; 0, a(3), a(2); a(3), 0, a(1)];
end
end

function out = BErs(b)
out = zeros(3, 6);
if isvector(b)
    out = [b(1), 0, 0, 1 / 2 * b(2), 0, 1 / 2 * b(3); 0, b(2), 0, 1 / 2 * b(1), 1 / 2 * b(3), 0; 0, 0, b(3), 0, 1 / 2 * b(2), 1 / 2 * b(1)];
end
end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -0.5 * (A(2, 3) + A(3, 2));
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -0.5 * (A(1, 3) + A(3, 1));
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 3) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 4) = -0.5 * A(3, 3);
D(4, 5) = 0.25 * (A(1, 3) + A(3, 1));
D(4, 6) = 0.25 * (A(2, 3) + A(3, 2));
D(5, 1) = -0.5 * (A(2, 3) + A(3, 2));
D(5, 4) = 0.25 * (A(1, 3) + A(3, 1));
D(5, 5) = -0.5 * A(1, 1);
D(5, 6) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 2) = -0.5 * (A(1, 3) + A(3, 1));
D(6, 4) = 0.25 * (A(2, 3) + A(3, 2));
D(6, 5) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 6) = -0.5 * A(2, 2);
end

function D = secDiffOperator2(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -A(2, 3);
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -A(1, 3);
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -A(1, 2);
D(4, 3) = -2 * A(1, 2);
D(4, 4) = -A(3, 3);
D(4, 5) = A(1, 3);
D(4, 6) = A(2, 3);
D(5, 1) = -2 * A(2, 3);
D(5, 4) = A(1, 3);
D(5, 5) = -A(1, 1);
D(5, 6) = A(1, 2);
D(6, 2) = -2 * A(1, 3);
D(6, 4) = A(2, 3);
D(6, 5) = A(1, 2);
D(6, 6) = -A(2, 2);
end