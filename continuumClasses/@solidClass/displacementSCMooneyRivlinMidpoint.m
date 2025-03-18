function [rData, kData, elementEnergy, array] = displacementSCMooneyRivlinMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
% 20.09.2016 Alexander Janz
% 19.08.2021 Marlon Franke: moofeKIT version

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;
numericalTangentObject = obj.numericalTangentObject;
% element degree of freedom tables and more
globalFullEdof = meshObject.globalFullEdof;
numberOfDOFs = size(globalFullEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
edof = meshObject.edof(e, :);
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

Re = rData{1};
Ke = kData{1, 1};
edN1 = dofs.edN1;
edR = qR(edof, 1:dimension)';
edN = qN(edof, 1:dimension)';
edN05 = 0.5 * (edN + edN1);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % Deformation gradient
    FN1 = edN1 * dN_X_I';
    FN = edN * dN_X_I';
    FN05 = edN05 * dN_X_I';
    % B-matrix (midpoint configuration)
    BN05 = BMatrix(dN_X_I, FN05);
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    CN = FN' * FN;
    CN05 = FN05' * FN05; %averaged
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1);
    GN05 = 0.5 * wedge(CN05, CN05); %mid-configuration
    % Third invariant
    cN1 = det(CN1);
    cN05 = det(CN05);
    % Derivative of the strain energy function
    DW_C = a * I;
    DW_G = b * I;
    if ~computePostData
        DW_c = -d / (2 * cN05) + c / 2 * (1 - 1 / sqrt(cN05));
 
        % Second Piola Kirchhoff stress tensor
        SN05 = 2 * (DW_C + wedge(DW_G, CN05) + DW_c * GN05);
        SN05_v = matrixToVoigt(SN05, 'stress');
%         SN05_v = [SN05(1, 1); SN05(2, 2); SN05(3, 3); SN05(1, 2); SN05(2, 3); SN05(1, 3)];

        % Residual
        Re = Re + BN05' * SN05_v * detJ * gaussWeight(k);

        % Tangent
%         D2W_c = d / (2 * cN05^2) + c / (2 * (cN05)^(3 / 2)); 
        D2W_c = d / (2 * cN05^2) + c / (4 * (cN05)^(3 / 2));

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

        D = CN05;
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
        Kmat3 = D2W_c * GN05_v * GN05_v';

        % Assembly of elasticity tensor
        ELA = (Kmat1 + Kmat2 + Kmat3);

        A1 = dN_X_I' * SN05 * dN_X_I * detJ * gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
%         Ke = Ke + BN05' * ELA * BN05 * detJ * gaussWeight(k) + 0.5 * MAT;
        Ke = Ke + 2 * BN05' * ELA * BN05 * detJ * gaussWeight(k) + 0.5 * MAT;
        % Strain energy function
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(cN1)) + c / 2 * (sqrt(cN1) - 1)^2) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        DW_cN1 = -d / (2 * cN1) + c / 2 * (1 - 1 / sqrt(cN1));
        SN1 = 2 * (DW_C + wedge(DW_G, CN1) + DW_cN1 * GN1);
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = Re;
    kData{1, 1} = Ke;
end
end