function [rData, kData, elementEnergy, array] = displacementMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Description:
%   material model: Mooney Rivlin in lambda and mu
%   with Psi_Vol = a*(log^2(J)+(J-1)^2) + d*log(J) and d = 2(a + 2b)
%   description in S and C
%   Integrator: Discrete Gradient
%
% 18.09.2023 Moritz Hille

analytic = true;

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% element degree of freedom tables and more
edof = meshObject.edof(e, :);
dimension = obj.dimension;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N_k_I;
dN = shapeFunctionObject.dN_xi_k_I;

% Voigt maps
% map.Voigt = [1, 5, 9, 4, 8, 7]';
% map.VoigtInv = [1, 4, 6; 4, 2, 5; 6, 5, 3];
% map.VoigtFull = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
% map.VoigtFullInv = [1, 4, 6; 7, 2, 5; 9, 8, 3];
% map.Isym = [eye(3), zeros(3); zeros(3), 2 * eye(3)];
% map.Isyminv = [eye(3), zeros(3); zeros(3), 1 / 2 * eye(3)];
% mapVoigtObject = obj.mapVoigtObject;

%initialize energy
% TODO

%% extract dofs for element
% nodal dofs
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
edR = qR(edof, 1:dimension)';
edN = qN(edof, 1:dimension)';
edN1 = qN1(edof, 1:dimension)';

%% material data
lambda = materialObject.lambda;
mu = materialObject.mu;

% parameter for Mooney Rivlin
a = mu / 2;
b = 2 * mu - a; % a+b=2*mu

%% initialize residual and tangent (global)
elementEnergy.internalEnergy = 0;
elementEnergy.deltaU = 0;

RG = rData{1}; % G = global = Residual for external (displacement) dofs
KGG = kData{1, 1};
dimR = numel(rData{1});

% assignment displacement dofs
correlation = zeros(dimR/dimension, 2); % auxiliary matrix for assignment of displacement dof for numerical tangent
for mm = 1:dimR / dimension
    correlation((3 * mm)-2, 1) = mm;
    correlation((3 * mm)-1, 1) = mm;
    correlation(3*mm, 1) = mm;
    correlation((3 * mm)-2, 2) = 1;
    correlation((3 * mm)-1, 2) = 2;
    correlation(3*mm, 2) = 3;
end

%% Run through all Gauss points
for k = 1:numberOfGausspoints

    KGG_GP = zeros(24);
    % Jacobian matrix and determinant
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dN, k, setupObject);

    % deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN = edN * dN_X_I';
    CxN1 = FxN1.' * FxN1;
    CxN = FxN.' * FxN;
    CxN1inv = CxN1 \ eye(dimension);
    CxNinv = CxN \ eye(dimension);

    % auxiliary variables
    JN1 = sqrt(det(CxN1));
    JN = sqrt(det(CxN));
    BN1 = BMatrix(dN_X_I, FxN1);

    if ~computePostData

        numTangEps = 1i * 1e-50; % increment for complex numerical tangent

        %% residual, energy and tangent

        % energy
        psiComN1 = internalEnergy(CxN1, JN1, lambda, a, b);
        psiComN = internalEnergy(CxN, JN, lambda, a, b);
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + psiComN1 * detJ * gaussWeight(k);
        elementEnergy.deltaU = elementEnergy.deltaU + (psiComN1 - psiComN) * gaussWeight(k);

        % Delta stress (threshhold MP vs DG)
        DeltaCx = CxN1 - CxN;
        DeltaCxV = matrixToVoigt(DeltaCx, 'stress');
        DeltaCxV3 = matrixToVoigt(DeltaCx, 'strain');
        normDeltaCx = norm(DeltaCx);

        if normDeltaCx < 100 * eps % case midpoint rule
            DG = false;

            % deformation gradient and auxiliary variables (midpoint)
            edN05 = 1 / 2 * (edN + edN1);
            FxN05 = edN05 * dN_X_I';
            CxN05 = FxN05.' * FxN05;
            CxN05inv = CxN05 \ eye(dimension);
            %             CxN05invV   = matrixToVoigt(CxN05inv, 'stress');
            JN05 = sqrt(det(CxN05));
            BN05 = BMatrix(dN_X_I, FxN05);

            % first derivatives of internal energy: e(C) = psiCom(C)
            DPsiComDCxN05 = firstDerivativeInternalEnergy(CxN05, CxN05inv, JN05, lambda, a, b);
            DDPsiComDCxN05DCxN05V = secondDerivativeInternalEnergy(CxN05inv, JN05, lambda, a, b);

            % stress
            SN05 = 2 * DPsiComDCxN05; % 2nd PK stress tensor midpoint rule
            SN05V = matrixToVoigt(SN05, 'stress'); % 2nd PK stress tensor in Voigt notation
            % residual
            RG = RG + BN05.' * SN05V * detJ * gaussWeight(k);
            % global tangent
            if analytic
                % analytic
                KGG_GP = kron(dN_X_I.'*DPsiComDCxN05*dN_X_I, eye(3)) + 2 * BN05.' * DDPsiComDCxN05DCxN05V * BN05;
            else
                % numeric
                for ii = 1:dimR
                    edN1Akt = edN1;
                    jj = correlation(ii, 2);
                    kk = correlation(ii, 1);
                    edN1Akt(jj, kk) = edN1Akt(jj, kk) + numTangEps;

                    %modified Cauchy Green
                    edN05 = 1 / 2 * (edN1Akt + edN);
                    FxN05 = edN05 * dN_X_I';
                    FxN1 = edN1Akt * dN_X_I';
                    FxN = edN * dN_X_I';
                    CxN05 = FxN05.' * FxN05;
                    CxN1 = FxN1.' * FxN1;
                    CxN = FxN.' * FxN;
                    CxN05inv = CxN05 \ eye(dimension);
                    CxN1inv = CxN1 \ eye(dimension);
                    CxNinv = CxN \ eye(dimension);
                    DeltaCx = CxN1 - CxN;

                    % B-matrix (midpoint)
                    BN05 = BMatrix(dN_X_I, FxN05);

                    % auxiliary variables
                    JN05 = sqrt(det(CxN05));
                    JN1 = sqrt(det(CxN1));
                    JN = sqrt(det(CxN));
                    % first derivatives of internal energy: psiCom(C)
                    DPsiComDCxN05 = firstDerivativeInternalEnergy(CxN05, CxN05inv, JN05, lambda, a, b);
                    % stress
                    SN05 = 2 * DPsiComDCxN05; % 2nd PK stress tensor midpoint rule
                    SN05V = matrixToVoigt(SN05, 'stress'); % 2nd PK stress tensor in Voigt notation
                    % numerical tangent
                    RGimag = BN05.' * SN05V;
                    KGG_GP(:, ii) = imag(RGimag) / imag(numTangEps);
                end
            end
        else % case discrete gradient
            DG = true;

            % deformation gradient and auxiliary variables (midpoint)
            edN05 = 1 / 2 * (edN + edN1);
            FxN05 = edN05 * dN_X_I';
            CxN05 = 1 / 2 * (CxN1 + CxN);
            CxN05inv = CxN05 \ eye(dimension);
            CxN05invV = matrixToVoigt(CxN05inv, 'stress');
            JN05 = sqrt(det(CxN05));
            BN05 = BMatrix(dN_X_I, FxN05);

            % first derivatives of internal energy: e(C) = psiCom(C)
            DPsiComDCxN1 = firstDerivativeInternalEnergy(CxN1, CxN1inv, JN1, lambda, a, b);
            DPsiComDCxN1V = matrixToVoigt(DPsiComDCxN1, 'stress');
            DPsiComDCxN05 = firstDerivativeInternalEnergy(CxN05, CxN05inv, JN05, lambda, a, b);
            DPsiComDCxN05V = matrixToVoigt(DPsiComDCxN05, 'stress');
            DDPsiComDCxN05DCxN05V = secondDerivativeInternalEnergy(CxN05inv, JN05, lambda, a, b);

            % stress
            DeltaCxsquared = dotdot(DeltaCx, DeltaCx);
            SAlgo = 2 * DPsiComDCxN05 + 2 * ((psiComN1 - psiComN - dotdot(DPsiComDCxN05, DeltaCx)) / DeltaCxsquared) * DeltaCx;
            SAlgoV = matrixToVoigt(SAlgo, 'stress'); % 2nd PK stress tensor in Voigt notation
            % residual
            RG = RG + BN05.' * SAlgoV * detJ * gaussWeight(k);
            % global tangent
            if analytic
                Isym = diag([1, 1, 1, 0.5, 0.5, 0.5]);
                CAlgo = 2 * DDPsiComDCxN05DCxN05V ...
                    +(1 / DeltaCxsquared) * (DeltaCxV * (4 * DPsiComDCxN1V - 4 * DPsiComDCxN05V - 2 * DDPsiComDCxN05DCxN05V * DeltaCxV3).') ...
                    +(2 / DeltaCxsquared) * (psiComN1 - psiComN - dotdot(DPsiComDCxN05, DeltaCx)) * (2 * Isym - 4 * (DeltaCxV * DeltaCxV.') * (1 / DeltaCxsquared));
                KGG_GP = kron(dN_X_I.'*DPsiComDCxN05*dN_X_I, eye(dimension)) + BN05.' * CAlgo * BN1;
            else
                for ii = 1:dimR
                    edN1Akt = edN1;
                    jj = correlation(ii, 2);
                    kk = correlation(ii, 1);
                    edN1Akt(jj, kk) = edN1Akt(jj, kk) + numTangEps;

                    %modified Cauchy Green
                    edN05 = 1 / 2 * (edN1Akt + edN);
                    FxN05 = edN05 * dN_X_I';
                    FxN1 = edN1Akt * dN_X_I';
                    FxN = edN * dN_X_I';
                    CxN05 = 1 / 2 * (CxN1 + CxN);
                    CxN1 = FxN1.' * FxN1;
                    CxN = FxN.' * FxN;
                    CxN05inv = CxN05 \ eye(dimension);
                    CxN1inv = CxN1 \ eye(dimension);
                    CxNinv = CxN \ eye(dimension);
                    DeltaCx = CxN1 - CxN;

                    % B-matrix (midpoint)
                    BN05 = BMatrix(dN_X_I, FxN05);

                    % auxiliary variables
                    JN05 = sqrt(det(CxN05));
                    JN1 = sqrt(det(CxN1));
                    JN = sqrt(det(CxN));

                    % energy
                    psiComN1 = internalEnergy(CxN1, JN1, lambda, a, b);
                    psiComN = internalEnergy(CxN, JN, lambda, a, b);

                    % first derivatives of internal energy: psiCom(C)
                    DPsiComDCxN05 = firstDerivativeInternalEnergy(CxN05, CxN05inv, JN05, lambda, a, b);
                    % stresses
                    SAlgo = 2 * DPsiComDCxN05 + 2 * ((psiComN1 - psiComN - dotdot(DPsiComDCxN05, DeltaCx)) / dotdot(DeltaCx, DeltaCx)) * DeltaCx;
                    SAlgoV = matrixToVoigt(SAlgo, 'stress'); % 2nd PK stress tensor in Voigt notation
                    %numerical tangent
                    RGimag = BN05.' * SAlgoV;
                    KGG_GP(:, ii) = imag(RGimag) / imag(numTangEps);
                end
            end
        end
        KGG = KGG + KGG_GP * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        epsilon = B * uN05(:);
        sigma_V = C * epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    DG;
    rData{1} = RG;
    kData{1, 1} = KGG;
end
end

%% auxiliary functions
function [psiVol] = internalEnergyVolumetricPart(J, lambda)
psiVol = (lambda / 2) * (log(J)^2 + (J - 1)^2);
end

function [energy] = internalEnergy(C, J, lambda, a, b)
traceCofC = 0.5 * (trace(C)^2 - trace(C*C));
energy = a * (trace(C) - 3 - 2 * log(J)) + b * (traceCofC - 3 - 4 * log(J)) + internalEnergyVolumetricPart(J, lambda);
end

function [f] = ftilde(J, lambda, a, b)
f = -a - 2 * b + (lambda / 2) * (log(J) + J^2 - J);
end

function [g] = g(J, lambda)
g = (lambda / 4) * (2 * J^2 - J + 1);
end

function [DPsi] = firstDerivativeInternalEnergy(C, Cinv, J, lambda, a, b)
DPsi = (a + b * trace(C)) * eye(3) - b * C + ftilde(J, lambda, a, b) * Cinv;
end

function [DDPsi] = secondDerivativeInternalEnergy(Cinv, J, lambda, a, b)
CinvV = matrixToVoigt(Cinv, 'stress');
unityV = [1; 1; 1; 0; 0; 0];
Isym = diag([1, 1, 1, 0.5, 0.5, 0.5]);
DDPsi = b * ((unityV * unityV.') - Isym) - ftilde(J, lambda, a, b) * sym_voigt(Cinv, Cinv, 'stress', 'stress') + g(J, lambda) * (CinvV * CinvV.');
end

function [dotdotAB] = dotdot(A, B)
dotdotAB = sum(sum(A.*B));
end

function [C_o_C_sym_voigt] = sym_voigt(A, B, stressStrainA, stressStrainB)
A11 = A(1, 1);
A12 = A(1, 2);
A13 = A(1, 3);
A21 = A(2, 1);
A22 = A(2, 2);
A23 = A(2, 3);
A31 = A(3, 1);
A32 = A(3, 2);
A33 = A(3, 3);

if strcmpi(stressStrainA, 'strain')
    A12 = 2 * A12;
    A21 = 2 * A21;
    A13 = 2 * A13;
    A31 = 2 * A31;
    A23 = 2 * A23;
    A32 = 2 * A32;
end

B11 = B(1, 1);
B12 = B(1, 2);
B13 = B(1, 3);
B21 = B(2, 1);
B22 = B(2, 2);
B23 = B(2, 3);
B31 = B(3, 1);
B32 = B(3, 2);
B33 = B(3, 3);

if strcmpi(stressStrainB, 'strain')
    B12 = 2 * B12;
    B21 = 2 * B21;
    B13 = 2 * B13;
    B31 = 2 * B31;
    B23 = 2 * B23;
    B32 = 2 * B32;
end

C_o_C_sym_voigt = [A11 * B11, A12 * B12, A13 * B13, 0.5 * A11 * B12 + 0.5 * A12 * B11, 0.5 * A12 * B13 + 0.5 * A13 * B12, 0.5 * A11 * B13 + 0.5 * A13 * B11; ...
    A21 * B21, A22 * B22, A23 * B23, 0.5 * A21 * B22 + 0.5 * A22 * B21, 0.5 * A22 * B23 + 0.5 * A23 * B22, 0.5 * A21 * B23 + 0.5 * A23 * B21; ...
    A31 * B31, A32 * B32, A33 * B33, 0.5 * A31 * B32 + 0.5 * A32 * B31, 0.5 * A32 * B33 + 0.5 * A33 * B32, 0.5 * A31 * B33 + 0.5 * A33 * B31; ...
    A11 * B21, A12 * B22, A13 * B23, 0.5 * A11 * B22 + 0.5 * A12 * B21, 0.5 * A12 * B23 + 0.5 * A13 * B22, 0.5 * A11 * B23 + 0.5 * A13 * B21; ...
    A21 * B31, A22 * B32, A23 * B33, 0.5 * A21 * B32 + 0.5 * A22 * B31, 0.5 * A22 * B33 + 0.5 * A23 * B32, 0.5 * A21 * B33 + 0.5 * A23 * B31; ...
    A11 * B31, A12 * B32, A13 * B33, 0.5 * A11 * B32 + 0.5 * A12 * B31, 0.5 * A12 * B33 + 0.5 * A13 * B32, 0.5 * A11 * B33 + 0.5 * A13 * B31];
end