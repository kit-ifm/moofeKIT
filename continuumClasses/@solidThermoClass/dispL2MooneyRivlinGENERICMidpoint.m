function [rData, kData, elementEnergy, array] = dispL2MooneyRivlinGENERICMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Description:
%   material model: Mooney Rivlin in a, b, c1, c2, d1 and d2 plus thermal evolution equation obtained from GENERIC
%   description in S, C and theta with mixed projection method for partial
%   derivatives of internal energy and entropy with respect to theta
%   Integrator: Midpoint Rule with L2-smoothing
%
% CREATOR(S)
% 04.02.2025 Moritz Hille, Tim Pleßke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

% element degree of freedom tables and more
edof = meshObject.edof(e, :);
dimension = obj.dimension;
DT = setupObject.timeStepSize;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N  = shapeFunctionObject.N_k_I;
dN = shapeFunctionObject.dN_xi_k_I;
M  = mixedFEObject.shapeFunctionObject.N_k_I;
dM = mixedFEObject.shapeFunctionObject.dN_xi_k_I;
numberOfInternalNodes = size(M, 2);

%% extract dofs for element
% nodal dofs
qR          = obj.qR;
qN          = obj.qN;
qN1         = dofs.edN1;
edR         = qR(edof, 1:dimension)';
edN         = qN(edof, 1:dimension)';
edN1        = qN1;
edN05       = 1 / 2 * (edN1 + edN);
vN05        = (edN1 - edN) / DT;
tempN1      = dofs.thetaN1;
tempN       = qN(edof,dimension+1).';
tempN05     = 1 / 2 * (tempN + tempN1);
edAlphaN    = obj.mixedFEObject.qN(e, :);
edAlphaN1   = dofs.edAlphaN1;
edAlphaN05  = 1 / 2 * (edAlphaN + edAlphaN1);

% mixed dofs
mixedEN     = edAlphaN(1:numberOfInternalNodes);
mixedEtaN   = edAlphaN(numberOfInternalNodes+1:2*numberOfInternalNodes);
mixedEN1    = edAlphaN1(1:numberOfInternalNodes);
mixedEtaN1  = edAlphaN1(numberOfInternalNodes+1:2*numberOfInternalNodes);
mixedEN05   = edAlphaN05(1:numberOfInternalNodes);
mixedEtaN05 = edAlphaN05(numberOfInternalNodes+1:2*numberOfInternalNodes);

%% material data
a       = materialObject.a;
b       = materialObject.b;
c1      = materialObject.c1;
c2      = materialObject.c2;
d1      = materialObject.d1;
d2      = materialObject.d2;

kappa   = materialObject.kappa;
beta    = materialObject.beta;
thetaR  = materialObject.thetaR;
k0      = materialObject.k0;

%% initialize energy, residual and tangent (global)

% initialize internal energy and difference of internal energy
elementEnergy.internalEnergy = 0;
elementEnergy.deltaU = 0;
elementEnergy.helmholtzEnergy = 0;
elementEnergy.entropy = 0;
elementEnergy.deltaS = 0;


% initialize sub-residual vector
RX   = rData{1};
RT   = rData{2};
RE   = rData{3};
REta = rData{4};

% initialize sub-tangent matrices (X = displacement part, T = temperature part)
KXX   = kData{1, 1}; KXT   = kData{1, 2}; KXE   = kData{1, 3}; KXEta   = kData{1, 4};
KTX   = kData{2, 1}; KTT   = kData{2, 2}; KTE   = kData{2, 3}; KTEta   = kData{2, 4};
KEX   = kData{3, 1}; KET   = kData{3, 2}; KEE   = kData{3, 3}; KEEta   = kData{3, 4};
KEtaX = kData{4, 1}; KEtaT = kData{4, 2}; KEtaE = kData{4, 3}; KEtaEta = kData{4, 4};

%% Run through all Gauss points
for k = 1:numberOfGausspoints

    % Jacobian matrix and determinant
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dN, k, setupObject);
    [~, ~, dM_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dM, k, setupObject);

    % deformation gradient, right cauchy green and velocity
    gradvN05    = vN05 * dN_X_I.';
    gradvN05V   = matrixToVoigtUnsymmetric(gradvN05);

    FxN1            = edN1  * dN_X_I.';
    FxN             = edN   * dN_X_I.';
    FxN05           = edN05 * dN_X_I.';
    FxNV            = matrixToVoigtUnsymmetric(FxN);
    FxN1V           = matrixToVoigtUnsymmetric(FxN1);
    FxN05V          = matrixToVoigtUnsymmetric(FxN05);
    FxN05inv        = FxN05 \ eye(dimension);
    FxN05invV       = matrixToVoigtUnsymmetric(FxN05inv);
    FxN05invtranspV = matrixToVoigtUnsymmetric(FxN05inv.');
    CxN1            = FxN1.'  * FxN1;
    CxN             = FxN.'   * FxN;
    CxN05           = FxN05.' * FxN05;

    CxNinv          = CxN \ eye(dimension);
    CxN1inv         = CxN1 \ eye(dimension);
    CxN05inv        = CxN05 \ eye(dimension);
    CxN05invV       = matrixToVoigt(CxN05inv, 'stress');

    MatTangTempV    = (CxN05invV * CxN05invV.') - 2*sym_voigt(CxN05inv, CxN05inv, 'stress', 'stress');

    % B-matrix (midpoint)
    BN05 = BMatrix(dN_X_I, FxN05);
    BUnsym = BMatrix(dN_X_I, 'mapType', 'unsymmetric');

    % auxiliary variables
    JN05 = sqrt(det(CxN05));
    JN1  = sqrt(det(CxN1 ));
    JN   = sqrt(det(CxN  ));

    % temperature
    thetaN1         = N(k, :) * tempN1.';
    thetaN          = N(k, :) * tempN.';
    thetaN05        = 1 / 2 * (thetaN + thetaN1);
    gradthetaN05    = tempN05 * dN_X_I.';

    % mixed variables
    EN       = M(k,:) * mixedEN.';
    EN1      = M(k,:) * mixedEN1.';
    EN05     = M(k,:) * mixedEN05.';
    EtaN     = M(k,:) * mixedEtaN.';
    EtaN1    = M(k,:) * mixedEtaN1.';
    EtaN05   = M(k,:) * mixedEtaN05.';
    gradEN05 = mixedEN05 * dM_X_I.';

    %% residual, energy and tangent

    % free energy (psi(C, G, J, theta) = psi0(theta) + psi1(C, G, J) + (theta - thetaR) * psi2(C))
    psi0N              = kappa * (thetaN   - thetaR - thetaN   * log(thetaN   / thetaR));
    psi0N1             = kappa * (thetaN1  - thetaR - thetaN1  * log(thetaN1  / thetaR));

    GxN                = JN^2*CxNinv;
    GxN1               = JN1^2*CxN1inv;

    psi1aN             = a * (trace(CxN) - 3) + b * (trace(GxN) - 3);
    psi1aN1            = a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3);

    psi1bN             = c1/2 * (JN-1)^2 - d1 * log(JN);
    psi1bN1            = c1/2 * (JN1-1)^2 - d1 * log(JN1);

    psi1N              = psi1aN + psi1bN;
    psi1N1             = psi1aN1 + psi1bN1;

    psi2N              = 3 * beta * (c2*(JN-1) - d2/(JN));
    psi2N1             = 3 * beta * (c2*(JN1-1) - d2/(JN1));

    % first derivatives free energy
    Dpsi0DthetaN       = -kappa * log(thetaN   / thetaR);
    Dpsi0DthetaN1      = -kappa * log(thetaN1  / thetaR);

    Dpsi1aDCxN1        = (a + b * trace(CxN1 )) * eye(dimension) -b * CxN1;
    Dpsi1aDCxN05       = (a + b * trace(CxN05)) * eye(dimension) -b * CxN05;

    Dpsi1bDCxN1        = (c1*(JN1 -1) - d1/(JN1 )) * 1/2 * JN1  * CxN1inv;
    Dpsi1bDCxN05       = (c1*(JN05-1) - d1/(JN05)) * 1/2 * JN05 * CxN05inv;

    Dpsi1DCxN1         = Dpsi1aDCxN1 + Dpsi1bDCxN1;
    Dpsi1DCxN05        = Dpsi1aDCxN05 + Dpsi1bDCxN05;

    Dpsi2DJN1          = 3 * beta * (c2 + d2 * JN1^(-2));
    Dpsi2DJN05         = 3 * beta * (c2 + d2 * JN05^(-2));

    Dpsi2DCxN1         = Dpsi2DJN1  * 1/2 * JN1  * CxN1inv;
    Dpsi2DCxN05        = Dpsi2DJN05 * 1/2 * JN05 * CxN05inv;

    % second derivatives free energy
    DDpsi2DJJN05       = - 6 * beta * d2 * JN05^(-3);

    % internal energy and entropy and their required derivatives
    internalEnergyN           = psi0N   - thetaN   * Dpsi0DthetaN   + psi1N   - thetaR * psi2N;
    internalEnergyN1          = psi0N1  - thetaN1  * Dpsi0DthetaN1  + psi1N1  - thetaR * psi2N1;
    entropyN                 = -Dpsi0DthetaN - psi2N;
    entropyN1                 = -Dpsi0DthetaN1 - psi2N1;

    DentropyDthetaN05         = kappa / thetaN05;

    DentropyDFxN05            = - Dpsi2DJN05 * JN05 * FxN05inv.';

    DentropyDCxN1             = - Dpsi2DJN1  * 1 / 2 * JN1  * CxN1inv;
    DentropyDCxN05            = - Dpsi2DJN05 * 1 / 2 * JN05 * CxN05inv;

    DentropyDFxN05V           = matrixToVoigtUnsymmetric(DentropyDFxN05);

    DDentropyDCxCxN05V        = -1 / 4 * DDpsi2DJJN05 * JN05^2 * (CxN05invV * CxN05invV.') - 1 / 2 * Dpsi2DJN05 * JN05 * (1 / 2 * (CxN05invV * CxN05invV.') - sym_voigt(CxN05inv, CxN05inv, 'stress', 'stress'));
    DDentropyDFxFxN05V        = -1 * ((DDpsi2DJJN05 * JN05^2 + Dpsi2DJN05 * JN05) * (FxN05invtranspV * FxN05invtranspV.') + Dpsi2DJN05 * (-JN05 * (FxN05invtranspV * FxN05invtranspV.') + DiffSecUnsym(FxN05)));
 
    DinternalEnergyDCxN1      = Dpsi1DCxN1  - thetaR * Dpsi2DCxN1;
    DinternalEnergyDCxN05     = Dpsi1DCxN05 - thetaR * Dpsi2DCxN05;

    unityV                    = [1;1;1;0;0;0];
    Isym                      = diag([1,1,1,0.5,0.5,0.5]);
    DDinternalEnergyDCxCxN05V = b * ((unityV*unityV.')- Isym) - 1/2 * (c1*JN05^2 - c1*JN05 - d1) * sym_voigt(CxN05inv, CxN05inv, 'stress', 'stress') + 1/4* JN05 * (2*c1*JN05 - c1)* (CxN05invV* CxN05invV.') + thetaR * DDentropyDCxCxN05V;

    % 2nd Piola Kirchhoff stress tensor and derivatives
    SN1            = 2 * (DinternalEnergyDCxN1  - thetaN1  * DentropyDCxN1);
    SN05           = 2 * (DinternalEnergyDCxN05 - thetaN05 * DentropyDCxN05);
    SN05V          = matrixToVoigt(SN05, 'stress');

    DSDCxN05V      = 2 * ( DDinternalEnergyDCxCxN05V - thetaN05 * DDentropyDCxCxN05V);
    DSDthetaN05    = -2 * DentropyDCxN05;
    DSDthetaN05V   = matrixToVoigt(DSDthetaN05, 'stress');

    if ~computePostData

        % heat flux vector
        QN05 = - k0 * JN05 * CxN05inv * gradthetaN05.';

        % residual
        RX   = RX   + BN05.' * SN05V * detJ * gaussWeight(k);
        RT   = RT   + (N(k,:).' * 1 / DT * (thetaN1 - thetaN) +  N(k,:).' * 1 / EtaN05 * (DentropyDFxN05V.' * gradvN05V) ...
            - ((1 / EN05 * dN_X_I.') - (N(k,:).' / EN05^2 * gradEN05)) * QN05) * detJ * gaussWeight(k);
        RE   = RE   + M(k,:).' * (kappa - EN05) * detJ * gaussWeight(k);
        REta = REta + M(k,:).' * (DentropyDthetaN05 - EtaN05) * detJ * gaussWeight(k);

        % analytical tangent
        KXX_GP     = kron(dN_X_I.' * 1 / 2 * SN05 * dN_X_I, eye(dimension)) + BN05.' * DSDCxN05V * BN05;
        KXT_GP     = BN05.' * 1 / 2 * DSDthetaN05V * N(k,:);
        %KXE_GP    = 0;
        %KXEta_GP  = 0;

        KTX_GP     = N(k,:).' / EtaN05 * (1 / 2 * gradvN05V.' * DDentropyDFxFxN05V + 1 / DT * DentropyDFxN05V.') * BUnsym ...
            + 1 / 2 * k0 * JN05 * (1 / EN05 * (BMatrixTemp(dN_X_I,gradthetaN05)).' - N(k,:).' / (EN05^2) * matrixToVoigt(1 / 2 * (gradEN05.' * gradthetaN05 + gradthetaN05.' * gradEN05),'strain').') * MatTangTempV * BN05;
        KTT_GP     = 1 / DT * N(k,:).' * N(k,:) + ((1 / EN05 * dN_X_I.') - (N(k,:).' / EN05^2 * gradEN05)) * k0 * JN05 * CxN05inv * 1 / 2 * dN_X_I;
        KTE_GP     = 1 / 2 * 1 / (EN05^2) * dN_X_I.' * QN05 * M(k,:) - N(k,:).' / (EN05^3) * gradEN05 * QN05 * M(k,:) + 1 / 2 * N(k,:).' / (EN05^2) * QN05.' * dM_X_I;
        KTEta_GP   = -1 / 2 * 1 / (EtaN05^2) * (DentropyDFxN05V.' * (gradvN05V)) * N(k,:).' * M(k,:);
        
        %KEX_GP    = 0;
        %KET_GP    = 0;
        KEE_GP     = -1 / 2 * M(k,:).' * M(k,:);
        %KEEta_GP  = 0;

        %KEtaX_GP  = 0;
        KEtaT_GP   = -1 / 2 * M(k,:).' * kappa / (thetaN05^2) * N(k,:);
        %KEtaE_GP  = 0;
        KEtaEta_GP = -1 / 2 * M(k,:).' * M(k,:);

        KXX        = KXX     + KXX_GP     * detJ * gaussWeight(k);
        KXT        = KXT     + KXT_GP     * detJ * gaussWeight(k);
        %KXE       = KXE     + KXE_GP     * detJ * gaussWeight(k);
        %KXEta     = KXEta   + KXEta_GP   * detJ * gaussWeight(k);

        KTX        = KTX     + KTX_GP     * detJ * gaussWeight(k);
        KTT        = KTT     + KTT_GP     * detJ * gaussWeight(k);
        KTE        = KTE     + KTE_GP     * detJ * gaussWeight(k);
        KTEta      = KTEta   + KTEta_GP   * detJ * gaussWeight(k);

        %KEX       = KEX     + KEX_GP     * detJ * gaussWeight(k);
        %KET       = KET     + KET_GP     * detJ * gaussWeight(k);
        KEE        = KEE     + KEE_GP     * detJ * gaussWeight(k);
        %KEEta     = KEEta   + KEEta_GP   * detJ * gaussWeight(k);

        %KEtaX     = KEtaX   + KEtaX_GP   * detJ * gaussWeight(k);
        KEtaT      = KEtaT   + KEtaT_GP   * detJ * gaussWeight(k);
        %KEtaE     = KEtaE   + KEtaE_GP   * detJ * gaussWeight(k);
        KEtaEta    = KEtaEta + KEtaEta_GP * detJ * gaussWeight(k);


        % energy and energy difference
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + internalEnergyN1 * detJ * gaussWeight(k);
        elementEnergy.deltaU = elementEnergy.deltaU + (internalEnergyN1 - internalEnergyN) * detJ * gaussWeight(k);
        % entropy and entropy difference
        elementEnergy.entropy = elementEnergy.entropy + entropyN1 * detJ * gaussWeight(k);
        elementEnergy.deltaS = elementEnergy.deltaS + (entropyN1 - entropyN) * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N, k, gaussWeight, detJ, stressTensor, setupObject, dimension);

    end
end

 if ~computePostData
    rData{1} = RX;
    rData{2} = RT;
    rData{3} = RE;
    rData{4} = REta;
    kData{1, 1} = KXX;   kData{1, 2} = KXT;   kData{1, 3} = KXE;   kData{1, 4} = KXEta;
    kData{2, 1} = KTX;   kData{2, 2} = KTT;   kData{2, 3} = KTE;   kData{2, 4} = KTEta;
    kData{3, 1} = KEX;   kData{3, 2} = KET;   kData{3, 3} = KEE;   kData{3, 4} = KEEta;
    kData{4, 1} = KEtaX; kData{4, 2} = KEtaT; kData{4, 3} = KEtaE; kData{4, 4} = KEtaEta;

 end
end

%% required functions

function Av = matrixToVoigtUnsymmetric(A)
%% Converts nonsymmetric matrix into voigt notation
Av = [A(1,1); A(2,2); A(3,3); A(1,2); A(2,3); A(1,3); A(2,1); A(3,2); A(3,1)];
end

function [C_o_C_sym_voigt] = sym_voigt(A,B, stressStrainA, stressStrainB)
A11 = A(1,1); A12 = A(1,2); A13 = A(1,3);
A21 = A(2,1); A22 = A(2,2); A23 = A(2,3);
A31 = A(3,1); A32 = A(3,2); A33 = A(3,3);

if strcmpi(stressStrainA,'strain')
    A12 = 2*A12;
    A21 = 2*A21;
    A13 = 2*A13;
    A31 = 2*A31;
    A23 = 2*A23;
    A32 = 2*A32;
end

B11 = B(1,1); B12 = B(1,2); B13 = B(1,3);
B21 = B(2,1); B22 = B(2,2); B23 = B(2,3);
B31 = B(3,1); B32 = B(3,2); B33 = B(3,3);

if strcmpi(stressStrainB,'strain')
    B12 = 2*B12;
    B21 = 2*B21;
    B13 = 2*B13;
    B31 = 2*B31;
    B23 = 2*B23;
    B32 = 2*B32;
end

C_o_C_sym_voigt = [ A11*B11   A12*B12   A13*B13  0.5*(A11*B12 + A12*B11)     0.5*(A12*B13 + A13*B12)     0.5*(A11*B13 + A13*B11);
    A21*B21   A22*B22   A23*B23  0.5*(A21*B22 + A22*B21)     0.5*(A22*B23 + A23*B22)     0.5*(A21*B23 + A23*B21);
    A31*B31   A32*B32   A33*B33  0.5*(A31*B32 + A32*B31)     0.5*(A32*B33 + A33*B32)     0.5*(A31*B33 + A33*B31);
    A11*B21   A12*B22   A13*B23  0.5*(A11*B22 + A12*B21)     0.5*(A12*B23 + A13*B22)     0.5*(A11*B23 + A13*B21);
    A21*B31   A22*B32   A23*B33  0.5*(A21*B32 + A22*B31)     0.5*(A22*B33 + A23*B32)     0.5*(A21*B33 + A23*B31);
    A11*B31   A12*B32   A13*B33  0.5*(A11*B32 + A12*B31)     0.5*(A12*B33 + A13*B32)     0.5*(A11*B33 + A13*B31)];
end

function DF = DiffSecUnsym(F)
F11 = F(1,1); F12 = F(1,2); F13 = F(1,3);
F21 = F(2,1); F22 = F(2,2); F23 = F(2,3);
F31 = F(3,1); F32 = F(3,2); F33 = F(3,3);

DF = [  0       F33     F22     0       -F32    0       0       -F23    0;
    F33     0       F11     0       0       -F31    0       0       -F13;
    F22     F11     0       -F21    0       0       -F12    0       0;
    0       0       -F21    0       F31     0       -F33    0       F23;
    -F32    0       0       F31     0       0       0       -F11    F12;
    0       -F31    0       0       0       0       F32     F21     -F22;
    0       0       -F12    -F33    0       F32     0       F13     0;
    -F23    0       0       0       -F11    F21     F13     0       0;
    0       -F13    0       F23     F12     -F22    0       0       0];
end

function [BQ] = BMatrixTemp(A, b)
BQ = zeros(6,8);
Ax = A(1,:);
Ay = A(2,:);
Az = A(3,:);
b1 = b(1);
b2 = b(2);
b3 = b(3);
BQ(1,:) = Ax * b1;
BQ(2,:) = Ay * b2;
BQ(3,:) = Az * b3;
BQ(4,:) = (Ax * b2 + Ay * b1);
BQ(5,:) = (Ay * b3 + Az * b2);
BQ(6,:) = (Ax * b3 + Az * b1);
end

% function A = voigtToMatrixUnsymmetric(Av)
% %% Converts nonsymmetric matrix into voigt notation
%     A = [Av(1) Av(4) Av(6);
%          Av(7) Av(2) Av(5);
%          Av(9) Av(8) Av(3)];
% end

% function [DC] = DiffSec(C)
%     C11 = C(1,1); C12 = C(1,2); C13 = C(1,3);
%     C21 = C(2,1); C22 = C(2,2); C23 = C(2,3);
%     C31 = C(3,1); C32 = C(3,2); C33 = C(3,3);
%
%     DC = [  0       C33     C22     0       -C32    0       0       -C23    0;
%             C33     0       C11     0       0       -C31    0       0       -C13;
%             C22     C11     0       -C21    0       0       -C12    0       0;
%             0       0       -C21    0       C31     0       -C33    0       C23;
%             -C32    0       0       C31     0       0       0       -C11    C12;
%             0       -C31    0       0       0       0       C32     C21     -C22;
%             0       0       -C12    -C33    0       C32     0       C13     0;
%             -C23    0       0       0       -C11    C21     C13     0       0;
%             0       -C13    0       C23     C12     -C22    0       0       0];
% end

% function [C_o_C_special_voigt] = special_voigt(A,B)
%     A11 = A(1,1); A12 = A(1,2); A13 = A(1,3);
%     A21 = A(2,1); A22 = A(2,2); A23 = A(2,3);
%     A31 = A(3,1); A32 = A(3,2); A33 = A(3,3);
%
%     B11 = B(1,1); B12 = B(1,2); B13 = B(1,3);
%     B21 = B(2,1); B22 = B(2,2); B23 = B(2,3);
%     B31 = B(3,1); B32 = B(3,2); B33 = B(3,3);
%
%     C_o_C_special_voigt = [ A11*B11     A12*B21     A13*B31     0.5*(A11*B21 + A12*B11)   0.5*(A12*B31 + A13*B21)   0.5*(A11*B31 + A13*B11)   0.5*(A11*B21 + A12*B11)   0.5*(A12*B31 + A13*B21)   0.5*(A11*B31 + A13*B11);
%                             A21*B12     A22*B22     A23*B32     0.5*(A21*B22 + A22*B12)   0.5*(A22*B32 + A23*B22)   0.5*(A21*B32 + A23*B12)   0.5*(A21*B22 + A22*B12)   0.5*(A22*B32 + A23*B22)   0.5*(A21*B32 + A23*B12);
%                             A31*B13     A32*B23     A33*B33     0.5*(A31*B23 + A32*B13)   0.5*(A32*B33 + A33*B23)   0.5*(A31*B33 + A33*B13)   0.5*(A31*B23 + A32*B13)   0.5*(A32*B33 + A33*B23)   0.5*(A31*B33 + A33*B13);
%                             A11*B12     A12*B22     A13*B32     0.5*(A11*B22 + A12*B12)   0.5*(A12*B32 + A13*B22)   0.5*(A11*B32 + A13*B12)   0.5*(A11*B22 + A12*B12)   0.5*(A12*B32 + A13*B22)   0.5*(A11*B32 + A13*B12);
%                             A21*B13     A22*B23     A23*B33     0.5*(A21*B23 + A22*B13)   0.5*(A22*B33 + A23*B23)   0.5*(A21*B33 + A23*B13)   0.5*(A21*B23 + A22*B13)   0.5*(A22*B33 + A23*B23)   0.5*(A21*B33 + A23*B13);
%                             A11*B13     A12*B23     A13*B33     0.5*(A11*B23 + A12*B13)   0.5*(A12*B33 + A13*B23)   0.5*(A11*B33 + A13*B13)   0.5*(A11*B23 + A12*B13)   0.5*(A12*B33 + A13*B23)   0.5*(A11*B33 + A13*B13);
%                             A21*B11     A22*B21     A23*B31     0.5*(A21*B21 + A22*B11)   0.5*(A22*B31 + A23*B21)   0.5*(A21*B31 + A23*B11)   0.5*(A21*B21 + A22*B11)   0.5*(A22*B31 + A23*B21)   0.5*(A21*B31 + A23*B11);
%                             A31*B12     A32*B22     A33*B32     0.5*(A31*B22 + A32*B12)   0.5*(A32*B32 + A33*B22)   0.5*(A31*B32 + A33*B12)   0.5*(A31*B22 + A32*B12)   0.5*(A32*B32 + A33*B22)   0.5*(A31*B32 + A33*B12);
%                             A31*B11     A32*B21     A33*B31     0.5*(A31*B21 + A32*B11)   0.5*(A32*B31 + A33*B21)   0.5*(A31*B31 + A33*B11)   0.5*(A31*B21 + A32*B11)   0.5*(A32*B31 + A33*B21)   0.5*(A31*B31 + A33*B11)];
% end

% function [C,CVoigt] = specialVoigtTensor(A,B)
%     C = zeros(3,3,3,3);
%     for I = 1:3
%         for J = 1:3
%             for K = 1:3
%                 for L = 1:3
%                     C(I,J,K,L) = 1/2*(A(I,K)*B(L,J)+A(I,L)*B(K,J));
%                 end
%             end
%         end
%     end
%     % TODO: build an operator for different Voigt matrix index notation instead of just one for the following fourth order Voigt tensor
%     CVoigt = [  C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,1,2) C(1,1,2,3) C(1,1,1,3) C(1,1,2,1) C(1,1,3,2) C(1,1,3,1);
%                 C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,1,2) C(2,2,2,3) C(2,2,1,3) C(2,2,2,1) C(2,2,3,2) C(2,2,3,1);
%                 C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,1,2) C(3,3,2,3) C(3,3,1,3) C(3,3,2,1) C(3,3,3,2) C(3,3,3,1);
%                 C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(1,2,1,2) C(1,2,2,3) C(1,2,1,3) C(1,2,2,1) C(1,2,3,2) C(1,2,3,1);
%                 C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(2,3,1,2) C(2,3,2,3) C(2,3,1,3) C(2,3,2,1) C(2,3,3,2) C(2,3,3,1);
%                 C(1,3,1,1) C(1,3,2,2) C(1,3,3,3) C(1,3,1,2) C(1,3,2,3) C(1,3,1,3) C(1,3,2,1) C(1,3,3,2) C(1,3,3,1);
%                 C(2,1,1,1) C(2,1,2,2) C(2,1,3,3) C(2,1,1,2) C(2,1,2,3) C(2,1,1,3) C(2,1,2,1) C(2,1,3,2) C(2,1,3,1);
%                 C(3,2,1,1) C(3,2,2,2) C(3,2,3,3) C(3,2,1,2) C(3,2,2,3) C(3,2,1,3) C(3,2,2,1) C(3,2,3,2) C(3,2,3,1);
%                 C(3,1,1,1) C(3,1,2,2) C(3,1,3,3) C(3,1,1,2) C(3,1,2,3) C(3,1,1,3) C(3,1,2,1) C(3,1,3,2) C(3,1,3,1)];
% end

% function [CVoigt,C] = A_obar_BVoigt(A,B)
%     C = zeros(3,3,3,3);
%     for I = 1:3
%         for J = 1:3
%             for K = 1:3
%                 for L = 1:3
%                     % C(I,J,K,L) = A(I,K)*B(L,J);
%                     % C(I,J,K,L) = A(I,K)*B(J,L);
%                     C(I,J,K,L) = -A(L,I)*B(J,K);
%                 end
%             end
%         end
%     end
%     CVoigt = [  C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,1,2) C(1,1,2,3) C(1,1,1,3) C(1,1,2,1) C(1,1,3,2) C(1,1,3,1);
%                 C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,1,2) C(2,2,2,3) C(2,2,1,3) C(2,2,2,1) C(2,2,3,2) C(2,2,3,1);
%                 C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,1,2) C(3,3,2,3) C(3,3,1,3) C(3,3,2,1) C(3,3,3,2) C(3,3,3,1);
%                 C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(1,2,1,2) C(1,2,2,3) C(1,2,1,3) C(1,2,2,1) C(1,2,3,2) C(1,2,3,1);
%                 C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(2,3,1,2) C(2,3,2,3) C(2,3,1,3) C(2,3,2,1) C(2,3,3,2) C(2,3,3,1);
%                 C(1,3,1,1) C(1,3,2,2) C(1,3,3,3) C(1,3,1,2) C(1,3,2,3) C(1,3,1,3) C(1,3,2,1) C(1,3,3,2) C(1,3,3,1);
%                 C(2,1,1,1) C(2,1,2,2) C(2,1,3,3) C(2,1,1,2) C(2,1,2,3) C(2,1,1,3) C(2,1,2,1) C(2,1,3,2) C(2,1,3,1);
%                 C(3,2,1,1) C(3,2,2,2) C(3,2,3,3) C(3,2,1,2) C(3,2,2,3) C(3,2,1,3) C(3,2,2,1) C(3,2,3,2) C(3,2,3,1);
%                 C(3,1,1,1) C(3,1,2,2) C(3,1,3,3) C(3,1,1,2) C(3,1,2,3) C(3,1,1,3) C(3,1,2,1) C(3,1,3,2) C(3,1,3,1)];
%     % sub2ind(size(C), [1 2], [1 2], [1 2], [1 2])
% end

% function C = OuterProduct(A, B)
%     C = reshape(A(:) * B(:).', [size(A), size(B)]);
% end

