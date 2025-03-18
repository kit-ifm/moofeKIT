function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinBatheDvorkinSkewSymbolicHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
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

% symbolic shape functions
syms xi eta;
syms xiBar etaBar;
syms w1 w2 w3 w4;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8;

load solFactors20231019.mat;
load testBsANS.mat;

x1e = ed(1, 1);
x2e = ed(1, 2);
x3e = ed(1, 3);
x4e = ed(1, 4);
y1e = ed(2, 1);
y2e = ed(2, 2);
y3e = ed(2, 3);
y4e = ed(2, 4);

% xNodes = [x1, x2, x3, x4];
% yNodes = [y1, y2, y3, y4];
% % 
N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);

x = N1*x1e+N2*x2e+N3*x3e+N4*x4e;
y = N1*y1e+N2*y2e+N3*y3e+N4*y4e;

J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);
J = [J11, J12; J21, J22];

JInv = inv(J.');

J0 = subs(J, {xi, eta}, {0, 0});
J011 = J0(1, 1);
J012 = J0(1, 2);
J021 = J0(2, 1);
J022 = J0(2, 2);
detJ0 = det(J0);

J0Inv = inv(J0.');
% % 
% % skew coordinates
% numberOfNodes = 4;
% nodesXi = elementNodesInLocalCoordinates(2, 'quadrilateral', 4);
% 
% [hourglassVector, ~] = computeHourglassVectorAndFunction(2, [xi; eta]);
% c = {0; 0};
% for ii = 1:numberOfNodes
%     c{1} = c{1} + 1 / numberOfNodes * xNodes(ii) * hourglassVector(ii);
%     c{2} = c{2} + 1 / numberOfNodes * yNodes(ii) * hourglassVector(ii);
% end
% H1 = xi * eta;
% pointsInSkewCoordinates = [xi; eta] + (J0 \ c) * H1;
% nodalPointsInSkewCoordinates = cell(2, 4);
% for ii = 1:numberOfNodes
%     nodalPointsInSkewCoordinates{1, ii} = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
%     nodalPointsInSkewCoordinates{2, ii} = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
% end
% 
% pointsInSkewCoordinates = [xiBar; etaBar];
% 
% % Metric shape functions
% xPoints = [nodalPointsInSkewCoordinates{1, :}].';
% yPoints = [nodalPointsInSkewCoordinates{2, :}].';
% P = [ones(4, 1), xPoints, yPoints, xPoints .* yPoints];
% p = [1, pointsInSkewCoordinates(1), pointsInSkewCoordinates(2), pointsInSkewCoordinates(1) .* pointsInSkewCoordinates(2)];
% M_k_I = p / P;
% M1 = M_k_I(1);
% M2 = M_k_I(2);
% M3 = M_k_I(3);
% M4 = M_k_I(4);
% 
% w = M1*w1+M2*w2+M3*w3+M4*w4;
% thetaX = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
% thetaY = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;
% wA = N1*w1+N2*w2+N3*w3+N4*w4;
% thetaXA = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
% thetaYA = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

% GammaXiBarO = subs(diff(w, xiBar)+J011*thetaX+J021*thetaY,[xiBar; etaBar] ,[xi; eta] + (J0 \ c) * H1);
% GammaEtaBarO = subs(diff(w, etaBar)+J012*thetaX+J022*thetaY,[xiBar; etaBar] ,[xi; eta] + (J0 \ c) * H1);
% [GammaXiBarC, GammaXiBarT] = coeffs(GammaXiBarO, [xi, eta]);
% [GammaEtaBarC, GammaEtaBarT] = coeffs(GammaEtaBarO, [xi, eta]);

% if length(GammaXiBarT) == 4
%     GammaXi = GammaXiBarC*(GammaXiBarT.*[alpha1, alpha2, alpha3, alpha4]).';
%     GammaEta = GammaEtaBarC*(GammaEtaBarT.*[alpha11, alpha12, alpha13, alpha14]).';
% else
%     GammaXi = GammaXiBarC*(GammaXiBarT.*[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7]).';
%     GammaEta = GammaEtaBarC*(GammaEtaBarT.*[alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17]).';
% end
% 

% 
% GammaXi = 1/4*(w2-w1+w3-w4)+1/8*((x2e-x1e)*(thetaX1+thetaX2)+(x3e-x4e)*(thetaX3+thetaX4))+1/8*((y2e-y1e)*(thetaY1+thetaY2)+(y3e-y4e)*(thetaY3+thetaY4));
% GammaEta = 1/4*(w4-w1+w3-w2)+1/8*((x4e-x1e)*(thetaX1+thetaX4)+(x3e-x2e)*(thetaX3+thetaX2))+1/8*((y4e-y1e)*(thetaY1+thetaY4)+(y3e-y2e)*(thetaY3+thetaY2));
% 
% GammaX = J0Inv(1, :) * [GammaXi; GammaEta];
% GammaY = J0Inv(2, :) * [GammaXi; GammaEta];
% 
% % Bending
% GammaXiBending = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
% GammaEtaBending = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
% 
% GammaXiBending2 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
% GammaEtaBending2 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
% 
% GammaXiBendingLinear = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);
% GammaEtaBendingLinear = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);
% 
% GammaXBending3 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
% GammaYBending3 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
% 
% GammaXiBendingLinear2 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);
% GammaEtaBendingLinear2 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);
% 
% GammaXiBendingLinear3 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1^2+y1^2+x1, x2*y2^2+y2^2+x2, x3*y3^2+y3^2+x3, x4*y4^2+y4^2+x4, -y1^2-1, -y2^2-1, -y3^2-1, -y4^2-1, -2*x1*y1-2*y1, -2*x2*y2-2*y2, -2*x3*y3-2*y3, -2*x4*y4-2*y4]);
% GammaEtaBendingLinear3 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1^2+y1^2, x2*y2^2+y2^2, x3*y3^2+y3^2, x4*y4^2+y4^2, -y1^2, -y2^2, -y3^2, -y4^2, -2*x1*y1-2*y1, -2*x2*y2-2*y2, -2*x3*y3-2*y3, -2*x4*y4-2*y4]);
% 
% % Shearing
% GammaXShearing = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
% GammaYShearing = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
% 
% GammaXShearing2 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);
% GammaYShearing2 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);
% 
% GammaXShearingLinear = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);
% GammaYShearingLinear = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);



% sol = solve([GammaXiBendingLinear; GammaEtaBendingLinear; GammaXiBendingLinear2; GammaEtaBendingLinear2; GammaXiBendingLinear3; GammaEtaBendingLinear3] == [0; 0; 0; 0; 0; 0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha11, alpha12, alpha13, alpha14, alpha15]);
% sol = vpasolve(eval(subs([GammaXiBending; GammaEtaBending; GammaXiBending2; GammaEtaBending2; GammaXShearing; GammaYShearing; GammaXShearing2; GammaYShearing2], [xi, eta], [gaussPointsInLocalCoordinates(1, 1), gaussPointsInLocalCoordinates(2, 1)])) == [0; 0; 0; 0; 1; 0; 0; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17]);

% sol = solve([GammaXiBending; GammaEtaBending; GammaXiBending2; GammaEtaBending2; GammaXShearing; GammaYShearing2] == [0; 0; 0; 0; 1; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17]);
% sol = solve([GammaXShearing; GammaYShearing; GammaXShearing2; GammaYShearing2] == [1; 0; 0; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17]);

% GammaXiCorrect = subs(GammaXi, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha11, sol.alpha12, sol.alpha13, sol.alpha14, sol.alpha15, sol.alpha16, sol.alpha17]);
% GammaEtaCorrect = subs(GammaEta, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha11, sol.alpha12, sol.alpha13, sol.alpha14, sol.alpha15, sol.alpha16, sol.alpha17]);
% 
% % if c{1} ~= 0 || c{2} ~= 0
% %     keyboard;
% % end
%



% [Ls1, Ls1T] = coeffs(GammaXi, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);
% [Ls2, Ls2T] = coeffs(GammaEta, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);
% 
% if length(Ls1T) == 8
%     LsANSSymbolic = sym(zeros(2, 12));
%     LsANSSymbolic(1, 1:3:end) = Ls1(1:2:end);
%     LsANSSymbolic(1, 2:3:end) = Ls1(2:2:end);
%     LsANSSymbolic(2, 1:3:end) = Ls2(1:2:end);
%     LsANSSymbolic(2, 3:3:end) = Ls2(2:2:end);
% else
%     LsANSSymbolic = sym(zeros(2, 12));
%     LsANSSymbolic(1, :) = Ls1;
%     LsANSSymbolic(2, :) = Ls2;
% end

KTest = zeros(2, 12);

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

    % shear component
    BsANSO = zeros(2, 3*numberOfNodes);
    BsANSO(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dNGamma_xi_k_I(1, 1, :) + (1 - gaussPoints(2, k)) * dNGamma_xi_k_I(1, 3, :));
    BsANSO(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dNGamma_xi_k_I(2, 4, :) + (1 - gaussPoints(1, k)) * dNGamma_xi_k_I(2, 2, :));
    BsANSO(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * NGamma_k_I(3, :));
    BsANSO(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * NGamma_k_I(2, :));
    BsANSO(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * NGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * NGamma_k_I(3, :));
    BsANSO(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * NGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * NGamma_k_I(2, :));
    BsANSSubs = eval(subs(subs(BsANS.', [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]), [xi, eta], [gaussPoints(1, k), gaussPoints(2, k)]));
    Bs = J.' \ BsANSSubs;

    %     LsANS = zeros(2, 3*numberOfNodes);
    %     LsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(2, k)) * J0(1, 1) * MGamma_k_I(1, :) + (1 - gaussPointsInSkewCoordinates(2, k)) * J0(1, 1) * MGamma_k_I(3, :));
    %     LsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(1, k)) * J0(1, 2) * MGamma_k_I(4, :) + (1 - gaussPointsInSkewCoordinates(1, k)) * J0(1, 2) * MGamma_k_I(2, :));
    %     LsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(2, k)) * J0(2, 1) * MGamma_k_I(1, :) + (1 - gaussPointsInSkewCoordinates(2, k)) * J0(2, 1) * MGamma_k_I(3, :));
    %     LsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPointsInSkewCoordinates(1, k)) * J0(2, 2) * MGamma_k_I(4, :) + (1 - gaussPointsInSkewCoordinates(1, k)) * J0(2, 2) * MGamma_k_I(2, :));
    %     Ls = J0.' \ LsANS;
    %     Ls(1, 1:3:end) = dM_X_I(1, :);
    %     Ls(2, 1:3:end) = dM_X_I(2, :);

%     LsANS = zeros(2, 3*numberOfNodes);
%     LsANS(1, 1:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * dMGamma_xi_A(1, :) + (1 - gaussPoints(2, k)) * dMGamma_xi_C(1, :));
%     LsANS(2, 1:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * dMGamma_xi_D(2, :) + (1 - gaussPoints(1, k)) * dMGamma_xi_B(2, :));
%     LsANS(1, 2:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(1, 1) * MGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(1, 1) * MGamma_k_I(3, :));
%     LsANS(2, 2:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(1, 2) * MGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(1, 2) * MGamma_k_I(2, :));
%     LsANS(1, 3:3:end) = 1 / 2 * ((1 + gaussPoints(2, k)) * JA(2, 1) * MGamma_k_I(1, :) + (1 - gaussPoints(2, k)) * JC(2, 1) * MGamma_k_I(3, :));
%     LsANS(2, 3:3:end) = 1 / 2 * ((1 + gaussPoints(1, k)) * JD(2, 2) * MGamma_k_I(4, :) + (1 - gaussPoints(1, k)) * JB(2, 2) * MGamma_k_I(2, :));
%     LsANS = eval(subs(LsANSSymbolic, [xi, eta], [gaussPointsInLocalCoordinates(1, k), gaussPointsInLocalCoordinates(2, k)]));

    LsANS = zeros(2, 12);
    LsANS(1, 2) = subs(sol_GammaXi_Constant.alpha5 + sol_GammaXi_Xi.alpha5 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha5 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha5 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 3) = subs(sol_GammaXi_Constant.alpha9 + sol_GammaXi_Xi.alpha9 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha9 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha9 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 5) = subs(sol_GammaXi_Constant.alpha6 + sol_GammaXi_Xi.alpha6 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha6 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha6 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 6) = subs(sol_GammaXi_Constant.alpha10 + sol_GammaXi_Xi.alpha10 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha10 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha10 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 8) = subs(sol_GammaXi_Constant.alpha7 + sol_GammaXi_Xi.alpha7 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha7 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha7 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 9) = subs(sol_GammaXi_Constant.alpha11 + sol_GammaXi_Xi.alpha11 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha11 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha11 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 11) = subs(sol_GammaXi_Constant.alpha8 + sol_GammaXi_Xi.alpha8 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha8 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha8 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(1, 12) = subs(sol_GammaXi_Constant.alpha12 + sol_GammaXi_Xi.alpha12 * gaussPoints(1, k) + sol_GammaXi_Eta.alpha12 * gaussPoints(2, k) + sol_GammaXi_XiEta.alpha12 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    
    LsANS(2, 2) = subs(sol_GammaEta_Constant.alpha5 + sol_GammaEta_Xi.alpha5 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha5 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha5 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 3) = subs(sol_GammaEta_Constant.alpha9 + sol_GammaEta_Xi.alpha9 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha9 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha9 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 5) = subs(sol_GammaEta_Constant.alpha6 + sol_GammaEta_Xi.alpha6 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha6 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha6 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 6) = subs(sol_GammaEta_Constant.alpha10 + sol_GammaEta_Xi.alpha10 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha10 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha10 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 8) = subs(sol_GammaEta_Constant.alpha7 + sol_GammaEta_Xi.alpha7 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha7 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha7 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 9) = subs(sol_GammaEta_Constant.alpha11 + sol_GammaEta_Xi.alpha11 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha11 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha11 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 11) = subs(sol_GammaEta_Constant.alpha8 + sol_GammaEta_Xi.alpha8 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha8 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha8 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);
    LsANS(2, 12) = subs(sol_GammaEta_Constant.alpha12 + sol_GammaEta_Xi.alpha12 * gaussPoints(1, k) + sol_GammaEta_Eta.alpha12 * gaussPoints(2, k) + sol_GammaEta_XiEta.alpha12 * gaussPoints(1, k) * gaussPoints(2, k), [x1, x2, x3, x4, y1, y2, y3, y4], [x1e, x2e, x3e, x4e, y1e, y2e, y3e, y4e]);

    Ls = eval(J0.' \ LsANS);
    Ls(1, 1:3:end) = dM_X_I(1, :);
    Ls(2, 1:3:end) = dM_X_I(2, :);


    LsO = zeros(2, 3*numberOfNodes);
    LsO(1, 1:3:end) = dM_X_I(1, :);
    LsO(2, 1:3:end) = dM_X_I(2, :);
    LsO(1, 2:3:end) = M_k_I(k, :);
    LsO(2, 3:3:end) = M_k_I(k, :);

    if ~computePostData
        % Tangent
        % bending component
        Kb = Kb + Bb' * Eb * Lb * detJ * gaussWeight(k);
        % shear component
        Ks = Ks + Bs' * Es * Ls * detJ * gaussWeight(k);

        KTest = KTest + (LsO-Ls)* detJ * gaussWeight(k);
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

%% RESIDUAL
RX = RX + KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end
