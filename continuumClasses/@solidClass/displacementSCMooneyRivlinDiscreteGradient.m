function [rData, kData, elementEnergy, array] = displacementSCMooneyRivlinDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
% 20.09.2016 Alexander Janz
% 19.08.2021 Marlon Franke: moofeKIT version

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfDOFs = size(globalFullEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
%   N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
qN = obj.qN;
% material data and voigt notation
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
I = eye(dimension);

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% initialize flagDiscreteGradient
if ~flagNumericalTangent
    initializeFlagDiscreteGradient(numericalTangentObject, shapeFunctionObject.numberOfGausspoints, 1)
end
flagDiscreteGradient = numericalTangentObject.flagDiscreteGradient;

Re = rData{1};
Ke = kData{1, 1};
edN1 = dofs.edN1;
edN = qN(edof(e, :), 1:dimension)';
edN05 = 0.5 * (edN + edN1);

J = qR(edof(e, :), 1:dimension)' * dNr';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, indx)') \ dNr(indx, :);
    % Deformation gradient
    FN1 = edN1 * dNx';
    FN = edN * dNx';
    FN05 = edN05 * dNx';
    % B-matrix (midpoint configuration)
    BN05 = BMatrix(dNx, FN05);
    BN05 = 2 * BN05;
    % B-matrix (current configuration)
    BN1 = BMatrix(dNx, FN1);
    BN1 = 2 * BN1;
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    CN = FN' * FN;
    CN05 = 0.5 * (CN + CN1); %averaged

    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1);
    GN = 0.5 * wedge(CN, CN);
    GN05 = 0.5 * (GN + GN1); %averaged
    GMid = 0.5 * wedge(CN05, CN05); %mid-configuration

    % Third invariant
    cN1 = det(CN1);
    cN = det(CN);
    cN05 = 0.5 * (cN + cN1);

    % Strain energy function
    elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2) * detJ * gaussWeight(k);

    % Derivative of the strain energy function
    DW_C = a * I;
    DW_G = b * I;
    deltac = cN1 - cN;
    normDeltac = abs(deltac);
    if ~flagNumericalTangent
        discreteGradientCondition_c = (normDeltac > 1e-10);
    else
        discreteGradientCondition_c = flagDiscreteGradient(k, 1);
    end
    if discreteGradientCondition_c
        %Greenspan formula '84
        DW_c = (-d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2 - (-d * log(sqrt(cN)) + c / 2 * (sqrt(cN) - 1)^2)) / (cN1 - cN);
        D2W_c = ((d / cN1 - (c * (cN1^(1 / 2) - 1)) / cN1^(1 / 2)) / (cN - cN1) + (c * (cN^(1 / 2) - 1)^2 - c * (cN1^(1 / 2) - 1)^2 - 2 * d * log(cN^(1 / 2)) + 2 * d * log(cN1^(1 / 2))) / (cN - cN1)^2);
    else
        DW_c = -d / (2 * cN05) + c / 2 * (1 - 1 / sqrt(cN05));
        D2W_c = d / (2 * cN05^2) + c / (2 * (cN05)^(3 / 2));
    end

    % Second Piola Kirchhoff stress tensor
    SN05 = 2 * (DW_C + wedge(DW_G, CN05) + DW_c * 1 / 3 * (wedge(CN05, CN05) + GN05));
    SN05_v = matrixToVoigt(SN05, 'stress');

    % Residual
    Re = Re + BN05' * 1 / 2 * SN05_v * detJ * gaussWeight(k);

    % Tangent

    % Derivative of wedge(Sigma_G,CN05)

    D = (DW_G);
    SecDiffOperator = [0, D(3, 3), D(2, 2), 0, -D(3, 2), 0; ...
        D(3, 3), 0, D(1, 1), 0, 0, -D(3, 1); ...
        D(2, 2), D(1, 1), 0, -D(2, 1), 0, 0; ...
        0, 0, -D(2, 1), -D(3, 3), D(3, 1), D(2, 3); ...
        -D(3, 2), 0, 0, D(3, 1), -D(1, 1), D(1, 2); ...
        0, -D(3, 1), 0, D(3, 2), D(2, 1), -D(2, 2)];
    SecDiffOperator(4:6, 4:6) = 0.5 * SecDiffOperator(4:6, 4:6);
    Kmat1 = SecDiffOperator;

    % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I

    D = (2 / 3 * (CN05 + 0.5 * CN1));
    SecDiffOperator = [0, D(3, 3), D(2, 2), 0, -D(3, 2), 0; ...
        D(3, 3), 0, D(1, 1), 0, 0, -D(3, 1); ...
        D(2, 2), D(1, 1), 0, -D(2, 1), 0, 0; ...
        0, 0, -D(2, 1), -D(3, 3), D(3, 1), D(2, 3); ...
        -D(3, 2), 0, 0, D(3, 1), -D(1, 1), D(1, 2); ...
        0, -D(3, 1), 0, D(3, 2), D(2, 1), -D(2, 2)];
    SecDiffOperator(4:6, 4:6) = 0.5 * SecDiffOperator(4:6, 4:6);

    Kmat2 = DW_c * SecDiffOperator;

    % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
    GN05_v = matrixToVoigt(GN05, 'stress');
    GMid_v = matrixToVoigt(GMid, 'stress');
    GN1_v = matrixToVoigt(GN1, 'stress');
    Kmat3 = D2W_c * ((GN05_v + 2 * GMid_v) * GN1_v') * 1 / 3;

    % Assembly of elasticity tensor
    ELA = (Kmat1 + Kmat2 + Kmat3);

    A1 = dNx' * SN05 * dNx * detJ * gaussWeight(k);
    MAT = zeros(numberOfDOFs);
    for g = 1:dimension
        MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
    end
    Ke = Ke + 0.5 * BN05' * ELA * BN1 * detJ * gaussWeight(k) + 0.5 * MAT;
end
if ~computePostData
    if ~flagNumericalTangent
        numericalTangentObject.flagDiscreteGradient = flagDiscreteGradient;
    end
    
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end