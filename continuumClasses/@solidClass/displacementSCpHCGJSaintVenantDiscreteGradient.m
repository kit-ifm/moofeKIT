function [rData, kData, elementEnergy, array] = displacementSCpHCGJSaintVenantDiscreteGradient(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
edof = meshObject.edof(e, :);
dimension = obj.dimension;

% aquire material data
lambda = materialObject.lambda;
mu = materialObject.mu;
DMat = [lambda + 2 * mu, lambda, lambda, 0, 0, 0; ...
    lambda, lambda + 2 * mu, lambda, 0, 0, 0; ...
    lambda, lambda, lambda + 2 * mu, 0, 0, 0; ...
    0, 0, 0, mu, 0, 0; ...
    0, 0, 0, 0, mu, 0; ...
    0, 0, 0, 0, 0, mu];

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;
edN05 = 0.5*(edN + edN1);

% compute nodal velocities
% edvN1 = obj.vN1(edof, 1:dimension).'; % does not work
DT = setupObject.timeStepSize;
% edvN1 = (edN1-edN)/DT;
edvN05 = (edN1-edN)/DT;

% initialize residual & tangent
RX = rData{1};
KXX = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % Deformation gradient
    FxN = edN * dN_X_I';
    FxN1 = edN1 * dN_X_I';
    FxN05 = edN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FxN05);
    % mixed/history variables
    FvN05 = edvN05 * dN_X_I';
    CvN05 = FxN05.'*FvN05 + FvN05.'*FxN05;
    % history variables
    CN = obj.historyN(e,k).C;
    CN1 = CN + DT*CvN05;
    % update history variables
    if ~flagNumericalTangent
        obj.historyN1(e,k).C = CN1;
    end
    % Derivative of the strain energy function
    EN_v = matrixToVoigt(1/2*(CN - eye(3)),'strain');
    EN1_v = matrixToVoigt(1/2*(CN1 - eye(3)),'strain');
    EAlgo_v = 0.5*(EN1_v+EN_v);
    % Second Piola Kirchhoff stress tensor
    SN1_v = DMat*EN1_v;
    SAlgo_v = DMat*EAlgo_v;
    if ~computePostData
        % Residual
        RX = RX + BN05.' * SAlgo_v * detJ * gaussWeight(k);
        % Strain energy
        W = 0.5 * EN1_v' * DMat * EN1_v;
        % W = a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3) + c / 2 * (JxN1 - 1)^2 - d * log(JxN1);
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
        PN1 = FxN1 * voigtToMatrix(SN1_v,'stress');
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
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
