function [rData, kData, elementEnergy, array] = displacementSCpHCGJMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
edN1 = dofs.edN1;

% compute nodal velocities
% edvN1 = obj.vN1(edof, 1:dimension).'; % does not work
DT = setupObject.timeStepSize;
edvN1 = (edN1-edN)/DT;

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
%     [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject);
    % Deformation gradient
    FxN1 = edN1 * dN_X_I';
    % B-matrix (current configuration)
    BN1 = BMatrix(dN_X_I, FxN1);
    % Right Cauchy-Green tensor
    CxN1 = FxN1.' * FxN1;
    % Cofactor
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    % Third invariant
    JxN1 = det(FxN1);
    % mixed/history variables
    FvN1 = edvN1 * dN_X_I';
    CvN1 = FxN1.'*FvN1 + FvN1.'*FxN1;
    GvN1 = wedge(CxN1,CvN1);
    HxN1 = 0.5 * wedge(FxN1, FxN1);
    CN = obj.historyN(e,k).C;
    GN = obj.historyN(e,k).G;
    JN = obj.historyN(e,k).J;  
    CN1 = CN + DT*CvN1;
    GN1 = GN + DT*GvN1;
    JN1 = JN + DT*HxN1(:).'*FvN1(:);
    % update history variables
    if ~flagNumericalTangent
        obj.historyN1(e,k).C = CN1;
        obj.historyN1(e,k).G = GN1;
        obj.historyN1(e,k).J = JN1;
    end
    % Derivative of the strain energy function
    dW_C = a * I;
    dW_G = b * I;
    % d_J_W = c*(JN1-1) - d*JN1^(-1);
    % d_J_d_J_W = c + d*JN1^(-2);
    dW_J = c*(JN1 - 1) - d*JN1^(-1);
    % Second Piola Kirchhoff stress tensor
    % SN1 = 2 * (dW_C + wedge(dW_G, CxN1) + 1/2*dW_J * JxN1*inv(CxN1));
    SN1 = 2 * (dW_C + wedge(dW_G, CxN1) + 1/2*dW_J*JxN1^(-1)*GxN1);
    SN1_v = matrixToVoigt(SN1, 'stress');
    if ~computePostData
        % Residual
        RX = RX + BN1.' * SN1_v * detJ * gaussWeight(k);
        % Strain energy
        W = a * (trace(CN1) - 3) + b * (trace(GN1) - 3) + c / 2 * (JN1 - 1)^2 - d * log(JN1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
%     disp(V)
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