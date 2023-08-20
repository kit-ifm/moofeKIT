function [rData, kData, elementEnergy, array] = displacementSCMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% displacementSCMooneyRivlinDiscreteGradient(obj,'PropertyName',PropertyValue)
%
% Description
% dispCascadeSC=: displacement based cascPade formulation
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 05.04.2017 M.S., A.J., M.F.

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
% mapVoigtObject = obj.mapVoigtObject;
% mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;

% general data
dimension = obj.dimension;
DT = setupObject.timeStepSize;

% dof data
numberOfDOFs = size(meshObject.globalFullEdof, 2);
mechanicalDOFs = 1:numberOfDOFs;
mechanicalDOFs(4:dimension+1:numberOfDOFs) = [];
numberOfMechanicalDOFs = numel(mechanicalDOFs);
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
edof = meshObject.edof;

% data for integration
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;

%% Create residual and tangent
% initialize residual, tangent & elementEnergy
RD = rData{1};
RT = rData{2};
KDD = kData{1, 1};
KDT = kData{1, 2};
KTD = kData{2, 1};
KTT = kData{2, 2};

elementEnergy.helmholtzEnergy = 0;
elementEnergy.internalEnergy = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 2)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

% extract dofs for element
edR = obj.qR(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
thetaN1 = dofs.thetaN1;
edN = obj.qN(edof(e, :), 1:dimension).';
thetaN = obj.qN(edof(e, :), dimension+1).';
edN05 = 1 / 2 * (edN + edN1);
thetaN05 = 1 / 2 * (thetaN + thetaN1);

JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % Temperature
    thetaNe = N_k_I(k, :) * thetaN.';
    thetaN1e = N_k_I(k, :) * thetaN1.';
    thetaN05e = 0.5 * (thetaN1e+thetaNe);

    % Deformation gradient
    FN = edN * dN_X_I';
    FN1 = edN1 * dN_X_I';
    FN05 = edN05 * dN_X_I';

    CN = FN.' * FN;
    CN1 = FN1.' * FN1;
    CN05 = FN05.' * FN05;
    CN05Average = 1 / 2 * (CN + CN1);

    GN = 0.5 * wedge(CN, CN);
    GN1 = 0.5 * wedge(CN1, CN1);
    GN05 = 0.5 * wedge(CN05, CN05);
    GN05Average = 0.5 * (GN1 + GN);

    cN = det(CN);
    cN1 = det(CN1);
    cN05 = det(CN05);

    % B-matrix
    BN1 = BMatrix(dN_X_I, FN1);

    % B-matrix (midpoint configuration)
    BN05 = BMatrix(dN_X_I, FN05);

    % Strain energy function
    PsiIsoN1 = a * (trace(CN1) - 3) + b * (trace(GN1) - 3);
    PsiVolN1 = -d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2;
    PsiThermoN1 = kappa * (thetaN1e-thetaR - thetaN1e * log(thetaN1e/thetaR));
    PsiCoupledN1 = -dimension * beta * (thetaN1e-thetaR) * (c * (sqrt(cN1) - 1) - d / sqrt(cN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    eThermN1 = kappa * (thetaN1e-thetaR);
    eCplN1 = 3 * beta * thetaR * (c * (sqrt(cN1) - 1) - d / sqrt(cN1));
    etotalN1 = PsiIsoN1 + PsiVolN1 + eThermN1 + eCplN1;
    elementEnergy.helmholtzEnergy = elementEnergy.helmholtzEnergy + PsiN1 * detJ * gaussWeight(k);
    elementEnergy.internalEnergy = elementEnergy.internalEnergy + etotalN1 * detJ * gaussWeight(k);

    % discrete Gradient conditions
    deltaC = CN1 - CN;
    normDeltaC = sqrt(innerProduct(deltaC, deltaC));
    deltaT = thetaN1e-thetaNe;
    normDeltaT = abs(deltaT);
    if ~flagNumericalTangent
        discreteGradientCondition_D = (normDeltaC > 10 * eps);
        discreteGradientCondition_T = (normDeltaT > 10 * eps);
    else
        discreteGradientCondition_D = flagDiscreteGradient(k, 1);
        discreteGradientCondition_T = flagDiscreteGradient(k, 2);
    end

    % Derivative of the strain energy function
    Du_C = a * eye(3);
    Du_G = b * eye(3);
    Du_theta = kappa;
    Deta_cN05 = 3 * beta * (c / (2 * sqrt(cN05)) + d / (2 * cN05^(3 / 2)));
    D2eta_c_cN05 = -3 * beta * (c / (4 * (cN05)^(3 / 2)) + 3 * d / (4 * (cN05)^(5 / 2)));
    if discreteGradientCondition_D
        flagDiscreteGradient(k, 1) = true;
        ucN1 = -d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2 + 3 * beta * thetaR * (c * (sqrt(cN1) - 1) - d / sqrt(cN1));
        ucN = -d * log(sqrt(cN)) + c / 2 * (sqrt(cN) - 1)^2 + 3 * beta * thetaR * (c * (sqrt(cN) - 1) - d / sqrt(cN));
        Du_c = (ucN1 - ucN) / (cN1 - cN);
        Du_c_c = (-d / (2 * cN1) + c / 2 * (1 - 1 / sqrt(cN1)) + 3 * beta * thetaR * (c / (2 * sqrt(cN1)) + d / (2 * (cN1)^(3 / 2)))) / (cN1 - cN) - (ucN1 - ucN) / (cN1 - cN)^2;
    else
        Du_c = -d / (2 * cN05) + c / 2 * (1 - 1 / sqrt(cN05)) + 3 * beta * thetaR * (c / 2 * 1 / sqrt(cN05) + d / 2 * cN05^(-3 / 2));
        Du_c_c = d / (2 * cN05^2) + c / (4 * cN05^(3 / 2)) - 3 * beta * thetaR * (c / (4 * cN05^(3 / 2)) + 3 * d / (4 * cN05^(5 / 2)));
    end
    if discreteGradientCondition_T
        flagDiscreteGradient(k, 2) = true;
        Deta_theta = (kappa * log(thetaN1e/thetaNe)) / (thetaN1e-thetaNe);
        D2eta_theta_theta = ((kappa / thetaN1e) * (thetaN1e-thetaNe) - kappa * log(thetaN1e/thetaNe)) / (thetaN1e-thetaNe)^2;
    else
        Deta_theta = kappa / thetaN05e;
        D2eta_theta_theta = -kappa / thetaN05e^2;
    end

    % create residual
    thetaAlgo = Du_theta * Deta_theta^(-1);
    GAlgo = 1 / 3 * (wedge(CN05Average, CN05Average) + GN05Average);
    if flagDiscreteGradient(k, 1)
        SAlgo = 2 * (Du_C + wedge(Du_G, CN05Average) + Du_c * GAlgo - thetaAlgo * Deta_cN05 * GN05);
    else
        SAlgo = 2 * (Du_C + wedge(Du_G, CN05) + (Du_c - thetaAlgo * Deta_cN05) * GN05);
    end

    RD = RD + BN05.' * matrixToVoigt(SAlgo, 'stress') * detJ * gaussWeight(k);

    QN05 = -k0 * (cN05^(-1) * GN05 * (dN_X_I * thetaN05.'));
    RT = RT + (N_k_I(k, :).' / DT * N_k_I(k, :) * (thetaN1 - thetaN).' + N_k_I(k, :).' * Deta_theta^(-1) * Deta_cN05 * innerProduct(GN05, 1/DT*(CN1 - CN)) - dN_X_I' * Du_theta^(-1) * QN05) * detJ * gaussWeight(k);

    % create tangent
    A1 = 1 / 2 * dN_X_I' * SAlgo * dN_X_I;
    KDDGeom = zeros(numberOfMechanicalDOFs);
    for g = 1:dimension
        KDDGeom(g:dimension:numberOfMechanicalDOFs, g:dimension:numberOfMechanicalDOFs) = A1;
    end
    if flagDiscreteGradient(k, 1)
        KDDMat1 = BN05' * 2 * secDiffOperator(Du_G) * BN1;
        KDDMat2 = BN05' * 4 * Du_c_c * matrixToVoigt(GAlgo, 'stress') * matrixToVoigt(GN1, 'stress')' * BN1;
        KDDMat3 = BN05' * 4 * Du_c * secDiffOperator(1/3*CN05Average+1/6*CN1) * BN1;
    else
        KDDMat1 = BN05' * 2 * secDiffOperator(Du_G) * BN05;
        KDDMat2 = BN05' * 2 * Du_c_c * matrixToVoigt(GN05, 'stress') * matrixToVoigt(GN05, 'stress')' * BN05;
        KDDMat3 = BN05' * 2 * Du_c * secDiffOperator(CN05) * BN05;
    end
    KDDMat4 = -BN05' * 2 * thetaAlgo * D2eta_c_cN05 * matrixToVoigt(GN05, 'stress') * matrixToVoigt(GN05, 'stress')' * BN05;
    KDDMat5 = -BN05' * 2 * thetaAlgo * Deta_cN05 * secDiffOperator(CN05) * BN05;
    KDD = KDD + (KDDGeom + KDDMat1 + KDDMat2 + KDDMat3 + KDDMat4 + KDDMat5) * detJ * gaussWeight(k);

    gamma = dN_X_I * thetaN05';
    gammaTransformed = [gamma(1), 0, 0, 1 / 2 * gamma(2), 0, 1 / 2 * gamma(3); 0, gamma(2), 0, 1 / 2 * gamma(1), 1 / 2 * gamma(3), 0; 0, 0, gamma(3), 0, 1 / 2 * gamma(2), 1 / 2 * gamma(1)];
    KTD1 = N_k_I(k, :).' * Deta_theta^(-1) * D2eta_c_cN05 * innerProduct(GN05, 1/DT*(CN1 - CN)) * matrixToVoigt(GN05, 'stress')' * BN05;
    KTD2 = N_k_I(k, :).' * Deta_theta^(-1) * Deta_cN05 * matrixToVoigt(wedge(CN05, 1/DT*(CN1 - CN)), 'stress')' * BN05;
    KTD3 = N_k_I(k, :).' * Deta_theta^(-1) * Deta_cN05 * matrixToVoigt(GN05, 'stress')' * 2 / DT * BN1;
    KTD4 = -dN_X_I' * Du_theta^(-1) * k0 * (cN05)^(-2) * GN05 * gamma * matrixToVoigt(GN05, 'stress')' * BN05;
    KTD5 = dN_X_I' * Du_theta^(-1) * k0 * (cN05)^(-1) * gammaTransformed * secDiffOperator2(CN05) * BN05;
    KTD = KTD + (KTD1 + KTD2 + KTD3 + KTD4 + KTD5) * detJ * gaussWeight(k);

    if flagDiscreteGradient(k, 2)
        KDT = KDT + kron(2*BN05'*Du_theta*Deta_theta^(-2)*D2eta_theta_theta*Deta_cN05*matrixToVoigt(GN05, 'stress'), N_k_I(k, :)) * detJ * gaussWeight(k);
        KTT = KTT + (N_k_I(k, :).' * 1 / DT * N_k_I(k, :) - N_k_I(k, :).' * Deta_theta^(-2) * D2eta_theta_theta * Deta_cN05 * innerProduct(GN05, 1/DT*(CN1 - CN)) * N_k_I(k, :) + 0.5 * dN_X_I' * Du_theta^(-1) * k0 * (cN05)^(-1) * GN05 * dN_X_I) * detJ * gaussWeight(k);
    else
        KDT = KDT + 0.5 * kron(2*BN05'*Du_theta*Deta_theta^(-2)*D2eta_theta_theta*Deta_cN05*matrixToVoigt(GN05, 'stress'), N_k_I(k, :)) * detJ * gaussWeight(k);
        KTT = KTT + (N_k_I(k, :).' * 1 / DT * N_k_I(k, :) - 0.5 * N_k_I(k, :).' * Deta_theta^(-2) * D2eta_theta_theta * Deta_cN05 * innerProduct(GN05, 1/DT*(CN1 - CN)) * N_k_I(k, :) + 0.5 * dN_X_I' * Du_theta^(-1) * k0 * (cN05)^(-1) * GN05 * dN_X_I) * detJ * gaussWeight(k);
    end
end
if ~computePostData
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
    rData{1} = RD;
    rData{2} = RT;
    kData{1, 1} = KDD;
    kData{1, 2} = KDT;
    kData{2, 1} = KTD;
    kData{2, 2} = KTT;
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