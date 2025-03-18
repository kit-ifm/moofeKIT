%% Quadrilateral
syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12 alpha13 alpha14;
syms beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8;

xNodes = [x1, x2, x3, x4];
yNodes = [y1, y2, y3, y4];

N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);

x = N1*x1+N2*x2+N3*x3+N4*x4;
y = N1*y1+N2*y2+N3*y3+N4*y4;

J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);
J = [J11, J12; J21, J22];
detJ = det(J);

JInv = inv(J.');

J0 = subs(J, {xi, eta}, {0, 0});
J011 = J0(1, 1);
J012 = J0(1, 2);
J021 = J0(2, 1);
J022 = J0(2, 2);
detJ0 = det(J);

J0Inv = inv(J0.');

% skew coordinates
numberOfNodes = 4;
nodesXi = elementNodesInLocalCoordinates(2, 'quadrilateral', 4);

[hourglassVector, ~] = computeHourglassVectorAndFunction(2, [xi; eta]);
c = {0; 0};
for ii = 1:numberOfNodes
    c{1} = c{1} + 1 / numberOfNodes * xNodes(ii) * hourglassVector(ii);
    c{2} = c{2} + 1 / numberOfNodes * yNodes(ii) * hourglassVector(ii);
end
H1 = xi * eta;
pointsInSkewCoordinates = [xi; eta] + (J0 \ c) * H1;
nodalPointsInSkewCoordinates = cell(2, 4);
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates{1, ii} = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates{2, ii} = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

pointsInSkewCoordinates = [xiBar; etaBar];

% Metric shape functions
xPoints = [nodalPointsInSkewCoordinates{1, :}].';
yPoints = [nodalPointsInSkewCoordinates{2, :}].';
P = [ones(4, 1), xPoints, yPoints, xPoints .* yPoints];
p = [1, pointsInSkewCoordinates(1), pointsInSkewCoordinates(2), pointsInSkewCoordinates(1) .* pointsInSkewCoordinates(2)];
M_k_I = p / P;
M1 = M_k_I(1);
M2 = M_k_I(2);
M3 = M_k_I(3);
M4 = M_k_I(4);

w = M1*w1+M2*w2+M3*w3+M4*w4;
thetaX = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
thetaY = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;

GammaXiBarO = subs(diff(w, xiBar)+J011*thetaX+J021*thetaY, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1);
GammaEtaBarO = subs(diff(w, etaBar)+J012*thetaX+J022*thetaY, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1);
[GammaXiBarC, GammaXiBarT] = coeffs(GammaXiBarO, [xi, eta]);
[GammaEtaBarC, GammaEtaBarT] = coeffs(GammaEtaBarO, [xi, eta]);
GammaXiBar = GammaXiBarC*(GammaXiBarT.*[alpha1, alpha2, alpha3, (alpha4+beta1+beta2), 0, 1, 1]).';
GammaEtaBar = GammaEtaBarC*(GammaEtaBarT.*[alpha8, alpha9, alpha10, (alpha11+beta3+beta4), 1, 0, 1]).';

GammaX = J0Inv(1, :) * [GammaXiBar; GammaEtaBar];
GammaY = J0Inv(2, :) * [GammaXiBar; GammaEtaBar];

% Bending
GammaXiBarBending = subs(GammaXiBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
GammaEtaBarBending = subs(GammaEtaBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);

GammaXiBarBendingRegular = subs(GammaXiBarBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBarBendingRegular = subs(GammaEtaBarBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBarBendingLinear = subs(GammaXiBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);
GammaEtaBarBendingLinear = subs(GammaEtaBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);

GammaXiBarBending2 = subs(GammaXiBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
GammaEtaBarBending2 = subs(GammaEtaBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);

GammaXiBarBending2Regular = subs(GammaXiBarBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBarBending2Regular = subs(GammaEtaBarBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBarBendingLinear2 = subs(GammaXiBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);
GammaEtaBarBendingLinear2 = subs(GammaEtaBar, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);


% Shearing
GammaXShearing = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXShearingRegular = subs(GammaXShearing, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaYShearingRegular = subs(GammaYShearing, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXShearing2 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing2 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXShearing2Regular = subs(GammaXShearing2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaYShearing2Regular = subs(GammaYShearing2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXShearingLinear = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearingLinear = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXShearingLinearRegular = subs(GammaXShearingLinear, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaYShearingLinearRegular = subs(GammaYShearingLinear, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

% solving
% sol = solve(subs([GammaXiBarBending; GammaXiBarBendingRegular; GammaXiBarBending2; GammaXiBarBending2Regular; GammaXShearing; GammaXShearingRegular; GammaXShearing2; GammaXShearing2Regular; GammaXShearingLinear; GammaXShearingLinearRegular; GammaEtaBending; GammaEtaBarBendingRegular; GammaEtaBarBending2; GammaEtaBarBending2Regular; GammaYShearing; GammaYShearingRegular; GammaYShearing2; GammaYShearing2Regular; GammaYShearingLinear; GammaYShearingLinearRegular], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]) == subs([0; 0; 0; 0; 1; 1; 0; 0; y; subs(y, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]); 0; 0; 0; 0; 0; 0; 1; 1; x; subs(x, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1])], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14]);
sol = solve(subs([GammaXiBarBending; GammaEtaBarBending; GammaXiBarBendingRegular; GammaEtaBarBendingRegular; GammaXiBarBending2; GammaEtaBarBending2; GammaXiBarBending2Regular; GammaEtaBarBending2Regular; GammaXShearing; GammaYShearing; GammaXShearingRegular; GammaYShearingRegular; GammaXShearing2; GammaYShearing2; GammaXShearing2Regular; GammaYShearing2Regular], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]) == subs([0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 1; 0; 0; 1; 0; 1], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha8, alpha9, alpha10, alpha11, beta1, beta2, beta3, beta4]);
% sol = solve([GammaXiBarBending; GammaEtaBarBending; GammaXiBarBendingRegular; GammaEtaBarBendingRegular; GammaXiBarBending2; GammaEtaBarBending2; GammaXiBarBending2Regular; GammaEtaBarBending2Regular; GammaXShearing; GammaYShearing; GammaXShearingRegular; GammaYShearingRegular; GammaXShearing2; GammaYShearing2; GammaXShearing2Regular; GammaYShearing2Regular] == subs([0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 1; 0; 0; 1; 0; 1], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14]);

% sol = solve([GammaXBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaYBending; GammaYBending2; GammaYShearing; GammaYShearing2;alpha2;alpha7] == [0; 0; 1; 0; 0; 0; 0; 1;0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
% solO = solve([GammaXBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaYBending; GammaYBending2; GammaYShearing; GammaYShearing2; GammaXBending3; GammaYBending3] == [0; 0; 1; 0; 0; 0; 0; 1; 0; 0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
GammaXiBarCorrect = subs(subs(GammaXiBar, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10, sol.alpha11, sol.alpha12, sol.alpha13, sol.alpha14]);
GammaEtaBarCorrect = subs(subs(GammaEtaBar, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10, sol.alpha11, sol.alpha12, sol.alpha13, sol.alpha14]), [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXCorrect = subs(GammaX, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]);
GammaYCorrect = subs(GammaY, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]);


