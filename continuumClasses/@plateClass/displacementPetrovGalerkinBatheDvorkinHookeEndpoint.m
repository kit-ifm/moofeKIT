function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinBatheDvorkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% DISPLACEMENTPETROVGALERKINBATHEDVORKINHOOKEENDPOINT Element routine of class plateClass.
%
% FORMULATION
% This is an Petrov-Galerkin Bathe-Dvorkin plate element.
% The standard lagrangian shape functions are employed for the virtual
% displacements while metric shape functions are used to approximate the
% displacement field. A Bathe-Dvorkin type ANS interpolation is applied to
% alleviate transverse shear locking.
% Only suitable for materially and geometrically linear simulations.
%
% CALL
% displacementPetrovGalerkinBatheDvorkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [1, 1] for residual data.
% kData: cell-array of size [1, 1] for tangent data.
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
% -
%
% SEE ALSO
% -
%
% CREATOR(S)
% Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

numberOfNodes = size(dN_xi_k_I, 3);
gaussPoints = shapeFunctionObject.gaussPoint;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% compute Jacobian matrices
J0 = ed * dN0_xi_I';
detJ0 = det(J0);
J011 = J0(1,1);
J012 = J0(1,2);
J021 = J0(2,1);
J022 = J0(2,2);

% compute material matrices
% bending component
Eb = zeros(3, 3);
Eb(1, 1) = 1;
Eb(1, 2) = nu;
Eb(2, 1) = nu;
Eb(2, 2) = 1;
Eb(3, 3) = (1 - nu) / 2;
Eb = E * h^3 / (12 * (1 - nu^2)) * Eb;
% shear component
Es = eye(2, 2);
Es = 5 / 6 * E / (2 * (1 + nu)) * h * Es;

% initialize energy
elementEnergy.strainEnergy = 0;

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

Kb = zeros(size(KXX));
Ks = zeros(size(KXX));

KsBubnov = zeros(size(KXX));

% compute Jacobian
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);

% shapeFunctions at collocation points
collocationPoints = [0, -1, 0, 1; 1, 0, -1, 0];
[NGamma_k_I, dNGamma_xi_k_I, ~] = computeLagrangeShapeFunction(2, 4, 4, collocationPoints);
JAllCollocationPoints = computeJacobianForAllGausspoints(ed, dNGamma_xi_k_I);
JA = extractJacobianForGausspoint(JAllCollocationPoints, 1, setupObject, dimension);
JB = extractJacobianForGausspoint(JAllCollocationPoints, 2, setupObject, dimension);
JC = extractJacobianForGausspoint(JAllCollocationPoints, 3, setupObject, dimension);
JD = extractJacobianForGausspoint(JAllCollocationPoints, 4, setupObject, dimension);

% Petrov: collocationPoints in skew
% Schubverzerrung mit metrischen Formfkten ausgewertet?!
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, ed, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, ed, J0, shapeFunctionObject);
collocationPointsInSkewCoordinates = computeSkewCoordinates(collocationPoints, ed, J0, shapeFunctionObject);
[MGamma_k_I, dMGamma_xiBar_k_I] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, collocationPointsInSkewCoordinates);

% compute additional shape functions
[M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);

% compute metric transverse shear strain
preFactorW = [nodesInSkewCoordinates(1, :); nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(1, :).^2; nodesInSkewCoordinates(2, :).^2; nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :).^2; zeros(1, 4)];
thetaXiVals = -[ones(1, 4); zeros(1, 4); nodesInSkewCoordinates(2, :); 2*nodesInSkewCoordinates(1, :); zeros(1, 4); 2*nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(2, :).^2; -nodesInSkewCoordinates(2, :)];
thetaEtaVals = -[zeros(1, 4); ones(1, 4); nodesInSkewCoordinates(1, :); zeros(1, 4); 2*nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(1, :).^2; 2*nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :); nodesInSkewCoordinates(1, :)];
thetaXVals = zeros(8, 4);
thetaYVals = zeros(8, 4);
for i=1:8
    thetaVals = inv(J0.')*[thetaXiVals(i, :); thetaEtaVals(i, :)];
    thetaXVals(i, :) = thetaVals(1, :);
    thetaYVals(i, :) = thetaVals(2, :);
end

shapeFunctionsGammaXi = inv([thetaXVals, thetaYVals])*([zeros(7, 4); gaussPointsInSkewCoordinates(2, :)]-preFactorW*dMr(1:2:end, :).');
shapeFunctionsGammaEta = inv([thetaXVals, thetaYVals])*([zeros(7, 4); -gaussPointsInSkewCoordinates(1, :)]-preFactorW*dMr(2:2:end, :).');


% Test: New Interpolation values
numberOfGausspoints1D = 2;
[gaussPoints1D, gaussWeights1D] = gaussPointsAndWeights(1, numberOfGausspoints1D, 'oneDimensional');
evalutionPoints1D = zeros(2, 4*numberOfGausspoints1D);
evalutionPoints1D(:, 1:numberOfGausspoints1D) = [gaussPoints1D; ones(1, numberOfGausspoints1D)];
evalutionPoints1D(:, numberOfGausspoints1D+1:2*numberOfGausspoints1D) = [-ones(1, numberOfGausspoints1D); gaussPoints1D];
evalutionPoints1D(:, 2*numberOfGausspoints1D+1:3*numberOfGausspoints1D) = [gaussPoints1D; -ones(1, numberOfGausspoints1D)];
evalutionPoints1D(:, 3*numberOfGausspoints1D+1:4*numberOfGausspoints1D) = [ones(1, numberOfGausspoints1D); gaussPoints1D];
evalutionPoints1DInSkewCoordinates = computeSkewCoordinates(evalutionPoints1D, ed, J0, shapeFunctionObject);
[MEvalutionPoints1D_k_I, ~] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, evalutionPoints1DInSkewCoordinates);

MGamma_A_I = zeros(1, 4);
MGamma_B_I = zeros(1, 4);
MGamma_C_I = zeros(1, 4);
MGamma_D_I = zeros(1, 4);
for k = 1:numberOfGausspoints1D
    MGamma_A_I = MGamma_A_I + MEvalutionPoints1D_k_I(k, :) * gaussWeights1D(k);
    MGamma_B_I = MGamma_B_I + MEvalutionPoints1D_k_I(k+numberOfGausspoints1D, :) * gaussWeights1D(k);
    MGamma_C_I = MGamma_C_I + MEvalutionPoints1D_k_I(k+2*numberOfGausspoints1D, :) * gaussWeights1D(k);
    MGamma_D_I = MGamma_D_I + MEvalutionPoints1D_k_I(k+3*numberOfGausspoints1D, :) * gaussWeights1D(k);
end

% Q Matrix
xVal = ed(1, :).';
yVal = ed(2, :).';
a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
h = 1/4*[1, -1, 1, -1].';
Q = zeros(8, 12);
Q(1, [2, 11]) = a1.'*xVal;
Q(1, [3, 12]) = a1.'*yVal;
Q(1, [5, 8]) = -a1.'*xVal;
Q(1, [6, 9]) = -a1.'*yVal;
Q(2, [2, 5]) = a2.'*xVal;
Q(2, [3, 6]) = a2.'*yVal;
Q(2, [8, 11]) = -a2.'*xVal;
Q(2, [9, 12]) = -a2.'*yVal;
Q(3, [5, 11]) = a1.'*xVal;
Q(3, [6, 12]) = a1.'*yVal;
Q(3, [2, 8]) = -a1.'*xVal;
Q(3, [3, 9]) = -a1.'*yVal;
Q(4, [5, 11]) = a2.'*xVal;
Q(4, [6, 12]) = a2.'*yVal;
Q(4, [2, 8]) = -a2.'*xVal;
Q(4, [3, 9]) = -a2.'*yVal;
Q(5, [5, 11]) = h.'*xVal; % TODO
Q(5, [6, 12]) = h.'*yVal;
Q(5, [2, 8]) = -h.'*xVal;
Q(5, [3, 9]) = -h.'*yVal;
Q(6, [5, 11]) = h.'*xVal; % TODO
Q(6, [6, 12]) = h.'*yVal;
Q(6, [2, 8]) = -h.'*xVal;
Q(6, [3, 9]) = -h.'*yVal;
Q(7, [5, 11]) = h.'*xVal; % TODO
Q(7, [6, 12]) = h.'*yVal;
Q(7, [2, 8]) = -h.'*xVal;
Q(7, [3, 9]) = -h.'*yVal;
Q(8, [5, 11]) = h.'*xVal; % TODO
Q(8, [6, 12]) = h.'*yVal;
Q(8, [2, 8]) = -h.'*xVal;
Q(8, [3, 9]) = -h.'*yVal;
Q = 1/4*Q;

C1 = 1/detJ0*(a2.'*yVal*h.'*xVal-a2.'*xVal*h.'*yVal);
C2 = 1/detJ0*(-a1.'*yVal*h.'*xVal+a1.'*xVal*h.'*yVal);

% test terms
termAS = 0;
termEASOrthStress = 0;
termEASOrthogonality = 0;
testOrtho = 0;

% test for BsTest
rightSide = zeros(4, 12);
leftTerm = zeros(4, 4);
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    BsU = zeros(2, 3*numberOfNodes);
    BsU(1, 1:3:end) = dN_X_I(1, :);
    BsU(2, 1:3:end) = dN_X_I(2, :);
    BsU(1, 2:3:end) = N_k_I(k, :);
    BsU(2, 3:3:end) = N_k_I(k, :);

    BsANS = zeros(2, 3*numberOfNodes);
    BsANS(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dNGamma_xi_k_I(1, 1, :) + (1 - gaussPoints(2, k)) * dNGamma_xi_k_I(1, 3, :));
    BsANS(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dNGamma_xi_k_I(2, 4, :) + (1 - gaussPoints(1, k)) * dNGamma_xi_k_I(2, 2, :));
    BsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * NGamma_k_I(3, :));
    BsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * NGamma_k_I(2, :));
    BsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * NGamma_k_I(3, :));
    BsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * NGamma_k_I(2, :));
    Bs = J.' \ BsANS;

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);
    xiBar = gaussPointsInSkewCoordinates(1, k);
    etaBar = gaussPointsInSkewCoordinates(2, k);

    S = [1, 0, etaBar, 0; 0, 1, 0, xiBar];
    EApprox = [xi^2*eta, 0, 1-eta^2, 0; 0, eta^2*xi, 0, 1-xi^2];
    EApprox = [xi, 0, xi*eta, 0; 0, eta, 0, xi*eta];
    EApprox = [xi*eta, 0, 1-eta^2, 0; 0, xi*eta, 0, 1-xi^2];

    rightSide = rightSide + (S.' * J0.' * BsU - S.' * J0.' * Bs) * detJ * gaussWeight(k);
    leftTerm = leftTerm + S.' * J0.' * inv(J.') * EApprox * detJ * gaussWeight(k);
end
% disp(rank(leftTerm));
BsTestPreFactors = inv(leftTerm)*rightSide;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    dMGamma_xi_A = (J0 \ JA).' * dMGamma_xiBar_k_I(1:2, :);
    dMGamma_xi_B = (J0 \ JB).' * dMGamma_xiBar_k_I(3:4, :);
    dMGamma_xi_C = (J0 \ JC).' * dMGamma_xiBar_k_I(5:6, :);
    dMGamma_xi_D = (J0 \ JD).' * dMGamma_xiBar_k_I(7:8, :);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dM_X_I = J0' \ dMr(indx, :);

    % B-Matrix
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);

    Lb = zeros(3, 3*numberOfNodes);
    Lb(1, 2:3:end) = dM_X_I(1, :);
    Lb(2, 3:3:end) = dM_X_I(2, :);
    Lb(3, 2:3:end) = dM_X_I(2, :);
    Lb(3, 3:3:end) = dM_X_I(1, :);

    % test solution for orthogonality constraint
    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);

    E = zeros(2, 6);
    E(1, 1) = gaussPoints(1, k);
    E(2, 2) = gaussPoints(2, k);
    E(1, 3) = gaussPoints(1, k)*gaussPoints(2, k)-1/3*C2*eta;
    E(2, 4) = gaussPoints(1, k)*gaussPoints(2, k)-1/3*C1*xi;
    E(1, 5) = 1/2*(3*gaussPoints(2, k).^2-1);
    E(2, 6) = 1/2*(3*gaussPoints(1, k).^2-1);
    E(1, 7) = gaussPoints(1, k)*gaussPoints(2, k).^2;
    E(2, 8) = gaussPoints(1, k).^2*gaussPoints(2, k);
    BsEps = detJ0/detJ*inv(J0.')*E*Q;
%     BsEps = detJ0/detJ*J0*E*Q;
    

    % shear component
    BsANS = zeros(2, 3*numberOfNodes);
    BsANS(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dNGamma_xi_k_I(1, 1, :) + (1 - gaussPoints(2, k)) * dNGamma_xi_k_I(1, 3, :));
    BsANS(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dNGamma_xi_k_I(2, 4, :) + (1 - gaussPoints(1, k)) * dNGamma_xi_k_I(2, 2, :));
    BsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * NGamma_k_I(3, :));
    BsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * NGamma_k_I(2, :));
    BsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * NGamma_k_I(3, :));
    BsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * NGamma_k_I(2, :));
    Bs = J.' \ BsANS;

    
    xiBar = gaussPointsInSkewCoordinates(1, k);
    etaBar = gaussPointsInSkewCoordinates(2, k);
    htx = h.'*xVal; % (x1/4 - x2/4 + x3/4 - x4/4);
    hty = h.'*yVal;

%     B = zeros(4, 12);
%     B(1, 1:3:end) = dNGamma_xi_k_I(1, 1, :);
%     B(2, 1:3:end) = dNGamma_xi_k_I(2, 2, :);
%     B(3, 1:3:end) = dNGamma_xi_k_I(1, 3, :);
%     B(4, 1:3:end) = dNGamma_xi_k_I(2, 4, :);
%     B(1, 2:3:end) = JA(1, 1) * NGamma_k_I(1, :);
%     B(2, 2:3:end) = JB(1, 2) * NGamma_k_I(2, :);
%     B(3, 2:3:end) = JC(1, 1) * NGamma_k_I(3, :);
%     B(4, 2:3:end) = JD(1, 2) * NGamma_k_I(4, :);
%     B(1, 3:3:end) = JA(2, 1) * NGamma_k_I(1, :);
%     B(2, 3:3:end) = JB(2, 2) * NGamma_k_I(2, :);
%     B(3, 3:3:end) = JC(2, 1) * NGamma_k_I(3, :);
%     B(4, 3:3:end) = JD(2, 2) * NGamma_k_I(4, :);
% 
%     BsNuevoBar = [eta/4 - 1/4, 2*(eta/2 - 1/2)*(xi/4 - 1/4)*(J011 - htx), 2*(eta/2 - 1/2)*(xi/4 - 1/4)*(J021 - hty),  1/4 - eta/4, -2*(eta/2 - 1/2)*(xi/4 + 1/4)*(J011 - htx), -2*(eta/2 - 1/2)*(xi/4 + 1/4)*(J021 - hty), eta/4 + 1/4, 2*(J011 + htx)*(eta/2 + 1/2)*(xi/4 + 1/4), 2*(J021 + hty)*(eta/2 + 1/2)*(xi/4 + 1/4), - eta/4 - 1/4, -2*(J011 + htx)*(eta/2 + 1/2)*(xi/4 - 1/4), -2*(J021 + hty)*(eta/2 + 1/2)*(xi/4 - 1/4)
%                  xi/4 - 1/4,   ((xi/2 - 1/2)*(J012 - htx)*(eta - 1))/2,   ((xi/2 - 1/2)*(J022 - hty)*(eta - 1))/2, - xi/4 - 1/4,   -((J012 + htx)*(xi/2 + 1/2)*(eta - 1))/2,   -((J022 + hty)*(xi/2 + 1/2)*(eta - 1))/2,  xi/4 + 1/4,   ((J012 + htx)*(xi/2 + 1/2)*(eta + 1))/2,   ((J022 + hty)*(xi/2 + 1/2)*(eta + 1))/2,    1/4 - xi/4,   -((xi/2 - 1/2)*(J012 - htx)*(eta + 1))/2,   -((xi/2 - 1/2)*(J022 - hty)*(eta + 1))/2];
% 
%     Bs = inv(J.')*BsNuevoBar;
    %     LsANS = zeros(2, 3*numberOfNodes);
    %     LsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(2, k)) * J0(1, 1) * MGamma_k_I(1, :) + (1 - gaussPointsInSkewCoordinates(2, k)) * J0(1, 1) * MGamma_k_I(3, :));
    %     LsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(1, k)) * J0(1, 2) * MGamma_k_I(4, :) + (1 - gaussPointsInSkewCoordinates(1, k)) * J0(1, 2) * MGamma_k_I(2, :));
    %     LsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(2, k)) * J0(2, 1) * MGamma_k_I(1, :) + (1 - gaussPointsInSkewCoordinates(2, k)) * J0(2, 1) * MGamma_k_I(3, :));
    %     LsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(1, k)) * J0(2, 2) * MGamma_k_I(4, :) + (1 - gaussPointsInSkewCoordinates(1, k)) * J0(2, 2) * MGamma_k_I(2, :));
    %     Ls = J0.' \ LsANS;
    %     Ls(1, 1:3:end) = dM_X_I(1, :);
    %     Ls(2, 1:3:end) = dM_X_I(2, :);

    LsANS = zeros(2, 3*numberOfNodes);
    LsANS(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dMGamma_xi_A(1, :) + (1 - gaussPoints(2, k)) * dMGamma_xi_C(1, :));
    LsANS(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dMGamma_xi_D(2, :) + (1 - gaussPoints(1, k)) * dMGamma_xi_B(2, :));
    LsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * MGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * MGamma_k_I(3, :));
    LsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * MGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * MGamma_k_I(2, :));
    LsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * MGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * MGamma_k_I(3, :));
    LsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * MGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * MGamma_k_I(2, :));
    Ls = J.' \ LsANS;

    A = zeros(4, 12);
    A(1, 1:3:end) = dMGamma_xi_A(1, :);
    A(2, 1:3:end) = dMGamma_xi_B(2, :);
    A(3, 1:3:end) = dMGamma_xi_C(1, :);
    A(4, 1:3:end) = dMGamma_xi_D(2, :);
    A(1, 2:3:end) = 1/2* JA(1, 1) * MGamma_A_I(1, :);
    A(2, 2:3:end) = 1/2* JB(1, 2) * MGamma_B_I(1, :);
    A(3, 2:3:end) = 1/2* JC(1, 1) * MGamma_C_I(1, :);
    A(4, 2:3:end) = 1/2* JD(1, 2) * MGamma_D_I(1, :);
    A(1, 3:3:end) = 1/2* JA(2, 1) * MGamma_A_I(1, :);
    A(2, 3:3:end) = 1/2* JB(2, 2) * MGamma_B_I(1, :);
    A(3, 3:3:end) = 1/2* JC(2, 1) * MGamma_C_I(1, :);
    A(4, 3:3:end) = 1/2* JD(2, 2) * MGamma_D_I(1, :);

    T = zeros(4, 4);
    T(1, 1) = JA(1, 1)*J0(2, 2)-J0(1,2)*JA(2,1);
    T(1, 2) = -JA(1, 1)*J0(2, 1)+J0(1,1)*JA(2,1);
    T(1, 3) = JA(1, 1)*J0(2, 2)-J0(1,2)*JA(2,1);
    T(2, 1) = JB(1, 2)*J0(2, 2)-J0(1,2)*JB(2,2);
    T(2, 2) = -JB(1, 2)*J0(2, 1)+J0(1,1)*JB(2,2);
    T(2, 4) = JB(1, 2)*J0(2, 1)-J0(1,1)*JB(2,2);
    T(3, 1) = JC(1, 1)*J0(2, 2)-J0(1,2)*JC(2,1);
    T(3, 2) = -JC(1, 1)*J0(2, 1)+J0(1,1)*JC(2,1);
    T(3, 3) = -JC(1, 1)*J0(2, 2)+J0(1,2)*JC(2,1);
    T(4, 1) = JD(1, 2)*J0(2, 2)-J0(1,2)*JD(2,2);
    T(4, 2) = -JD(1, 2)*J0(2, 1)+J0(1,1)*JD(2,2);
    T(4, 4) = -JD(1, 2)*J0(2, 1)+J0(1,1)*JD(2,2);
    T = 1/detJ0*T;

    LsANS = [1, 0, gaussPointsInSkewCoordinates(2, k), 0; 0, 1, 0, gaussPointsInSkewCoordinates(1, k)]*inv(T)*A;
    Ls = J0.' \ LsANS;

%     LsDSG = [-(((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)),                       ((((J011 - htx)*(C1^2*C2 + 2*C1^2 + C1*C2^2 + C2^2 - C2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J011 - htx)*(- C1^2*C2 - C1*C2^2 + C2^2 + C2))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J011 + htx)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1*C2^2 + 2*C1*C2 + C2^2 + C2)*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), ((((J021 - hty)*(C1^2*C2 + 2*C1^2 + C1*C2^2 + C2^2 - C2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J021 - hty)*(- C1^2*C2 - C1*C2^2 + C2^2 + C2))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J021 + hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1*C2^2 + 2*C1*C2 + C2^2 + C2)*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)),                   - ((((J011 - htx)*(- C1^2*C2 + C1*C2^2 - C2^2 + C2))/(6*(C1^2 + C2^2 - 1)) - ((J011 - htx)*(- C1^2*C2 + 2*C1^2 + C1*C2^2 + C2^2 + C2 - 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J011 + htx)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1*C2^2 + 2*C1*C2 - C2^2 + C2)*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), - ((((J021 - hty)*(- C1^2*C2 + C1*C2^2 - C2^2 + C2))/(6*(C1^2 + C2^2 - 1)) - ((J021 - hty)*(- C1^2*C2 + 2*C1^2 + C1*C2^2 + C2^2 + C2 - 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J021 + hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1*C2^2 + 2*C1*C2 - C2^2 + C2)*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)),                       ((((J011 + htx)*(- C1^2*C2 + 2*C1^2 - C1*C2^2 + C2^2 + C2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J011 + htx)*(C1^2*C2 + C1*C2^2 + C2^2 - C2))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J011 - htx)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1*C2^2 - 2*C1*C2 - C2^2 + C2)*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), ((((J021 + hty)*(- C1^2*C2 + 2*C1^2 - C1*C2^2 + C2^2 + C2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J021 + hty)*(C1^2*C2 + C1*C2^2 + C2^2 - C2))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J021 - hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1*C2^2 - 2*C1*C2 - C2^2 + C2)*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), -(((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)),                     ((((J011 + htx)*(- C1^2*C2 + C1*C2^2 + C2^2 + C2))/(6*(C1^2 + C2^2 - 1)) - ((J011 + htx)*(- C1^2*C2 - 2*C1^2 + C1*C2^2 - C2^2 + C2 + 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J011 - htx)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1*C2^2 - 2*C1*C2 + C2^2 + C2)*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2), ((((J021 + hty)*(- C1^2*C2 + C1*C2^2 + C2^2 + C2))/(6*(C1^2 + C2^2 - 1)) - ((J021 + hty)*(- C1^2*C2 - 2*C1^2 + C1*C2^2 - C2^2 + C2 + 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1 - etaBar + C1*etaBar + C2*etaBar + C2^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J021 - hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1*C2^2 - 2*C1*C2 + C2^2 + C2)*(etaBar - C1 + C1*etaBar - C2*etaBar + C2^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2); ...
%    -(((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)), (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*((J021*(- C1^2*C2 + C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)) - (hty*(2*C1^2 + 2*C2^2 - 2))/(6*(C1^2 + C2^2 - 1)) + (J021*(C1^2*C2 + C1^2 + C1*C2^2 - C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + (J021*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1^2 + C1*C2^2 + 2*C1*C2 + C1)*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),       ((((J022 - hty)*(C1^2*C2 + C1^2 + C1*C2^2 - C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J022 - hty)*(- C1^2*C2 + C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J022 + hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 + C1^2 + C1*C2^2 + 2*C1*C2 + C1)*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),   -(((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)), (J021*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(- C1^2*C2 + C1^2 + C1*C2^2 - 2*C1*C2 + C1)*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2) - (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*((hty*(2*C1^2 + 2*C2^2 - 2))/(6*(C1^2 + C2^2 - 1)) - (J021*(C1^2*C2 + C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)) + (J021*(C1^2*C2 - C1^2 - C1*C2^2 + C1 - 2*C2^2 + 2))/(2*(C1^2 + C2^2 - 1)))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)),           ((((J022 + hty)*(C1^2*C2 + C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)) - ((J022 + hty)*(C1^2*C2 - C1^2 - C1*C2^2 + C1 - 2*C2^2 + 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) + ((J022 - hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(- C1^2*C2 + C1^2 + C1*C2^2 - 2*C1*C2 + C1)*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),    (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)), (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*((hty*(2*C1^2 + 2*C2^2 - 2))/(6*(C1^2 + C2^2 - 1)) + (J021*(- C1^2*C2 + C1^2 - C1*C2^2 + C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)) + (J021*(C1^2*C2 + C1^2 + C1*C2^2 - C1))/(6*(C1^2 + C2^2 - 1)))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - (J021*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1^2 + C1*C2^2 - 2*C1*C2 + C1)*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),       ((((J022 + hty)*(- C1^2*C2 + C1^2 - C1*C2^2 + C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)) + ((J022 + hty)*(C1^2*C2 + C1^2 + C1*C2^2 - C1))/(6*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J022 - hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(C1^2*C2 - C1^2 + C1*C2^2 - 2*C1*C2 + C1)*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),     (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)), (((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*((hty*(2*C1^2 + 2*C2^2 - 2))/(6*(C1^2 + C2^2 - 1)) + (J021*(C1^2*C2 + C1^2 - C1*C2^2 + C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)) - (J021*(C1^2*C2 - C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - (J021*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(- C1^2*C2 - C1^2 + C1*C2^2 + 2*C1*C2 + C1)*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2),       - ((((J022 - hty)*(C1^2*C2 - C1^2 - C1*C2^2 + C1))/(6*(C1^2 + C2^2 - 1)) - ((J022 - hty)*(C1^2*C2 + C1^2 - C1*C2^2 + C1 + 2*C2^2 - 2))/(2*(C1^2 + C2^2 - 1)))*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(xiBar - C2 - C1*xiBar + C2*xiBar + C1^2 + C1*C2 - 1))/(4*(C1^2 + C2^2 - 1)) - ((J022 + hty)*((J011*J022)/(J011*J022 - J012*J021) - (J012*J021)/(J011*J022 - J012*J021))*(- C1^2*C2 - C1^2 + C1*C2^2 + 2*C1*C2 + C1)*(C2 - xiBar + C1*xiBar + C2*xiBar + C1^2 - C1*C2 - 1))/(12*(C1^2 + C2^2 - 1)^2)];
%     Ls = J0.' \ LsDSG;

%     BsANS = [1, 0, gaussPointsInSkewCoordinates(2, k), 0; 0, 1, 0, gaussPointsInSkewCoordinates(1, k)]*inv(T)*B;
%     Bs = J0.' \ BsANS;

    BsU = zeros(2, 3*numberOfNodes);
    BsU(1, 1:3:end) = dN_X_I(1, :);
    BsU(2, 1:3:end) = dN_X_I(2, :);
    BsU(1, 2:3:end) = N_k_I(k, :);
    BsU(2, 3:3:end) = N_k_I(k, :);

    LsU = zeros(2, 3*numberOfNodes);
    LsU(1, 1:3:end) = dM_X_I(1, :);
    LsU(2, 1:3:end) = dM_X_I(2, :);
    LsU(1, 2:3:end) = M_k_I(k, :);
    LsU(2, 3:3:end) = M_k_I(k, :);

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);
    EApprox = [xi, 0, xi*eta, 0; 0, eta, 0, xi*eta];
    BsTest = inv(J.') * EApprox * BsTestPreFactors;
    

    if ~computePostData
        % Tangent
        % bending component
        Kb = Kb + Bb' * Eb * Lb * detJ * gaussWeight(k);
        % shear component
        Ks = Ks + ((Bs).' * Es * Ls)* detJ * gaussWeight(k);

        KsBubnov = KsBubnov + ((Bs).' * Es * Bs)* detJ * gaussWeight(k);

        % test terms
        S = J0*[1, 0, gaussPointsInSkewCoordinates(2, k), 0; 0, 1, 0, gaussPointsInSkewCoordinates(1, k)];
        termAS = termAS + (((BsU) - (Bs + BsTest)).' * S) * detJ * gaussWeight(k);
        termEASOrthStress = termEASOrthStress + ((BsU+BsEps).' * Es * (Ls-LsU))* detJ * gaussWeight(k);
        termEASOrthogonality = termEASOrthogonality + (BsEps.' * Es * Ls)* detJ * gaussWeight(k);
        testOrtho = testOrtho + E.' * adjoint(J0)*Ls;
    else
        % stress at gausspoint
        kappa = Lb * qN1(:);
        gamma = Ls * qN1(:);

        m = Eb * kappa;
        q = Es * gamma;
        stressTensor.Cauchy = [m(1), m(3), q(1); m(3), m(2), q(2); q(1), q(2), 0];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

KXX = KXX + Kb + Ks;
% disp(max(max(abs(termAS))));

% test Ks
wTest = nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :).^2;
thetaXiTest = -2*nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :).^2;
thetaEtaTest = -2*nodesInSkewCoordinates(2, :).*nodesInSkewCoordinates(1, :).^2;
thetaXYTest = inv(J0.')*[thetaXiTest; thetaEtaTest];
thetaXTest = thetaXYTest(1, :);
thetaYTest = thetaXYTest(2, :);
qTest = [wTest, thetaXTest, thetaYTest];
qTestFormat = [qTest(1:4:end),qTest(2:4:end), qTest(3:4:end),qTest(4:4:end)];
testKs = Ks*qTestFormat.';

wTest = nodesInSkewCoordinates(1, :).^2;
thetaXiTest = -2*nodesInSkewCoordinates(1, :);
thetaEtaTest = zeros(1, 4);
thetaXYTest = inv(J0.')*[thetaXiTest; thetaEtaTest];
thetaXTest = thetaXYTest(1, :);
thetaYTest = thetaXYTest(2, :);
qTest = [wTest, thetaXTest, thetaYTest];
qTestFormat = [qTest(1:4:end),qTest(2:4:end), qTest(3:4:end),qTest(4:4:end)];
testKs = Ks*qTestFormat.';

wTest = nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :);
thetaXiTest = -nodesInSkewCoordinates(2, :);
thetaEtaTest = -nodesInSkewCoordinates(1, :);
thetaXYTest = inv(J0.')*[thetaXiTest; thetaEtaTest];
thetaXTest = thetaXYTest(1, :);
thetaYTest = thetaXYTest(2, :);
qTest = [wTest, thetaXTest, thetaYTest];
qTestFormat = [qTest(1:4:end),qTest(2:4:end), qTest(3:4:end),qTest(4:4:end)];
testKs = Ks*qTestFormat.';

wTest = nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :);
thetaXiTest = -2*nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :);
thetaEtaTest = -ones(1, 4);
thetaXYTest = inv(J0.')*[thetaXiTest; thetaEtaTest];
thetaXTest = thetaXYTest(1, :);
thetaYTest = thetaXYTest(2, :);
qTest = [wTest, thetaXTest, thetaYTest];
qTestFormat = [qTest(1:4:end),qTest(2:4:end), qTest(3:4:end),qTest(4:4:end)];
testKs = Ks*qTestFormat.';

%% RESIDUAL
RX = RX + KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end
