function [rData, kData, elementEnergy, array] = mixedSCpHCGJMooneyRivlinVol2Midpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCMOONEYRIVLINENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic Mooney-Rivlin
% ('MooneyRivlin') model (attention: nonlinear geometric and material/
% stress-strain relation).
% The formulation is based on a polyconvexity inspired strain-energy
% density function.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementSCMooneyRivlinEndpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2018.01.013
%
% SEE ALSO
% displacementSCMooneyRivlinMidpoint,
% displacementSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

DT = setupObject.timeStepSize;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);
numberOfXDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalNodes = size(M_k_I, 2);

dimension = obj.dimension;

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5*(edN + edN1);

% compute nodal velocities
% edvN1 = obj.vN1(edof, 1:dimension).'; % does not work
edvN05 = (edN1-edN)/DT;

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

edAlphaN = mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedCNv = edAlphaN(1:6*numberOfInternalNodes);
extractedGNv = edAlphaN(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN = edAlphaN(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes);
extractedGN1v = edAlphaN1(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes);

% initialize residual
RX = rData{1};
RDotCv = rData{2};
RDotGv = rData{3};
RDotJ = rData{4};

% initialize tangent
% KXX = kData{1, 1};
% KCC = kData{2, 2};
% KCG = kData{2, 3};
% KGC = kData{3, 2};
% KJJ = kData{4, 4};

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
%     [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % compute the values of the variables at the current Gauss point
    CNv = reshape(extractedCNv, 6, []) * M_k_I(k, :)';
    CN = voigtToMatrix(CNv, 'stress');
    CN1v = reshape(extractedCN1v, 6, []) * M_k_I(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    CN05 = 0.5 * (CN + CN1);

    GNv = reshape(extractedGNv, 6, []) * M_k_I(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M_k_I(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    GN05 = 0.5 * (GN + GN1);

    JN = extractedJN.' * M_k_I(k, :)';
    JN1 = extractedJN1.' * M_k_I(k, :)';
    JN05 = 0.5 * (JN + JN1);

    % Deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FxN05);

    % Right Cauchy-Green tensor
    CxN1 = FxN1.' * FxN1;
    CxN05 = FxN05.' * FxN05;
    % Cofactor
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    GxN05 = 0.5 * wedge(CxN05, CxN05);
    % Determinant
    JxN1 = det(FxN1);
    JxN05 = det(FxN05);
    % Derivative of the strain energy function
    dW_C = a * I;
    dW_G = b * I;
    % dW_J = c*(JN05 - 1) - d*JN05^(-1);
    dW_J = c*(JN05 - 1) - d*JN05;
    % d_J_d_J_W = c + d*JN1^(-2);
    % Second Piola Kirchhoff stress tensor
    SN1 = 2 * (dW_C + wedge(dW_G, CxN1) + 1/2*dW_J*JxN1^(-1)*GxN1);
    SN05 = 2 * (dW_C + wedge(dW_G, CxN05) + 1/2*dW_J*JxN05^(-1)*GxN05);
    FvN05 = edvN05 * dN_X_I';
    CvN05 = FxN05.'*FvN05 + FvN05.'*FxN05;
    GvN05 = wedge(CxN05,CvN05);
    HxN05 = 0.5 * wedge(FxN05, FxN05);
    SN05_v = matrixToVoigt(SN05, 'stress');
    if ~computePostData
        % Residual
        RX = RX + BN05.' * SN05_v * detJ * gaussWeight(k);
        RDotCv = RDotCv + kron(M_k_I(k, :)',matrixToVoigt((CN1 - CN) - DT*CvN05,'stress')) * detJ * gaussWeight(k);
        RDotGv = RDotGv + kron(M_k_I(k, :)',matrixToVoigt((GN1 - GN) - DT*GvN05,'stress')) * detJ * gaussWeight(k);
        RDotJ = RDotJ + M_k_I(k, :)' * ((JN1 - JN) - DT*HxN05(:).'*FvN05(:)) * detJ * gaussWeight(k);
        % Strain energy
        W = a * (trace(CN1) - 3) + b * (trace(GN1) - 3) + c / 2 * (JN1 - 1)^2 - d * log(JN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
        
        % % Tangent
        % Sigma_I_I = c/2*JxN1^(-2) + d*JxN1^(-3);    
        % % Derivative of wedge(Sigma_G,CN05)
        % if flagSymmetric
        %     Kmat1 = secDiffOperator(Sigma_G);
        % else
        %     Kmat1 = b*[  0     1     1     0     0     0     0     0     0;
        %                  1     0     1     0     0     0     0     0     0;
        %                  1     1     0     0     0     0     0     0     0;
        %                  0     0     0     0     0     0    -1     0     0;
        %                  0     0     0     0     0     0     0    -1     0;
        %                  0     0     0     0     0     0     0     0    -1;
        %                  0     0     0    -1     0     0     0     0     0;
        %                  0     0     0     0    -1     0     0     0     0;
        %                  0     0     0     0     0    -1     0     0     0];
        % end
        % % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I
        % if flagSymmetric
        %     Kmat2 = Sigma_I * secDiffOperator(CN1);
        % else
        %     Kmat2 = Sigma_I*[   0     CN1(3,3)     CN1(2,2)     0     -CN1(3,2)     0     0     -CN1(2,3)     0;
        %                          CN1(3,3)     0     CN1(1,1)     0     0     -CN1(3,1)     0     0     -CN1(1,3);
        %                          CN1(2,2)     CN1(1,1)     0     -CN1(2,1)     0     0     -CN1(1,2)     0     0;
        %                          0     0     -CN1(2,1)     0     CN1(3,1)     0    -CN1(3,3)     0     CN1(2,3);
        %                          -CN1(3,2)     0     0     CN1(3,1)     0     0     0    -CN1(1,1)     CN1(1,2);
        %                          0     -CN1(3,1)     0     0     0     0     CN1(3,2)     CN1(2,1)    -CN1(2,2);
        %                          0     0     -CN1(1,2)    -CN1(3,3)     0     CN1(3,2)     0     CN1(1,3)     0;
        %                          -CN1(2,3)     0     0     0    -CN1(1,1)     CN1(2,1)     CN1(1,3)     0     0;
        %                          0     -CN1(1,3)     0     CN1(2,3)     CN1(1,2)    -CN1(2,2)     0     0     0];
        % end
        % % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
        % if flagSymmetric
        %     GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3)];
        % else
        %     GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3); GN1(2, 1); GN1(3, 2); GN1(3, 1)];
        % end
        % Kmat3 = Sigma_I_I * 1/2 * JN1^(-1) * (GN1_v * GN1_v');
        % % Assembly of elasticity tensor
        % ELA = 4 * (Kmat1 + Kmat2 + Kmat3);
        % MAT = BN1' * ELA * BN1 * detJ * gaussWeight(k);
        % GEO = kron(dN_X_I' * SN1 * dN_X_I * detJ * gaussWeight(k),eye(3));
        % KXX = KXX + MAT + GEO;
    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    % pass residual
    rData{1} = RX;
    rData{2} = RDotCv;
    rData{3} = RDotGv;
    rData{4} = RDotJ;

    % pass tangent
    % kData{1, 1} = KXX;
    % kData{2, 2} = KCC;
    % kData{2, 3} = KCG;
    % kData{3, 2} = KGC;
    % kData{4, 4} = KJJ;
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