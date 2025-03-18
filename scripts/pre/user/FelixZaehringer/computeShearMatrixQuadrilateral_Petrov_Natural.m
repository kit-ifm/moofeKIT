%% Quadrilateral
syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12 alpha13 alpha14 alpha15 alpha16 alpha17 alpha18 alpha19 alpha20;

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

% pointsInSkewCoordinates = [xiBar; etaBar];

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

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;
[GammaXiC, GammaXiT] = coeffs(GammaXiO, [xi, eta]);
[GammaEtaC, GammaEtaT] = coeffs(GammaEtaO, [xi, eta]);
GammaXi = GammaXiC*(GammaXiT.*[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10]).';
GammaEta = GammaEtaC*(GammaEtaT.*[alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20]).';

GammaX = JInv(1, :) * [GammaXi; GammaEta];
GammaY = JInv(2, :) * [GammaXi; GammaEta];

% Bending
GammaXiBending = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
GammaEtaBending = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);

GammaXiBendingRegular = subs(GammaXiBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBendingRegular = subs(GammaEtaBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBending2 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
GammaEtaBending2 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);

GammaXiBending2Regular = subs(GammaXiBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBending2Regular = subs(GammaEtaBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBendingLinear = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);
GammaEtaBendingLinear = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);

GammaXiBendingLinearRegular = subs(GammaXiBendingLinear, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBendingLinearRegular = subs(GammaEtaBendingLinear, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

% GammaXBending3 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);
% GammaYBending3 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2, 1/2*x2^2, 1/2*x3^2, 1/2*x4^2, -x1, -x2, -x3, -x4, 0, 0, 0, 0]);

GammaXiBendingLinear2 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);
GammaEtaBendingLinear2 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1*y1^2, 1/2*x2*y2^2, 1/2*x3*y3^2, 1/2*x4*y4^2, -1/2*y1^2, -1/2*y2^2, -1/2*y3^2, -1/2*y4^2, -x1*y1, -x2*y2, -x3*y3, -x4*y4]);

GammaXiBendingLinear3 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1^2, x2*y2^2, x3*y3^2, x4*y4^2, -y1^2, -y2^2, -y3^2, -y4^2, -2*x1*y1, -2*x2*y2, -2*x3*y3, -2*x4*y4]);
GammaEtaBendingLinear3 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1^2, x2*y2^2, x3*y3^2, x4*y4^2, -y1^2, -y2^2, -y3^2, -y4^2, -2*x1*y1, -2*x2*y2, -2*x3*y3, -2*x4*y4]);

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


GammaXiShearing = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaEtaShearing = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXiShearingLinear = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaEtaShearingLinear = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1, x2*y2, x3*y3, x4*y4, 0, 0, 0, 0, 0, 0, 0, 0]);

% solving
sol = solve(subs([GammaXiBending; GammaXiBendingRegular; GammaXiBending2; GammaXiBending2Regular; GammaXShearing; GammaXShearingRegular; GammaXShearing2; GammaXShearing2Regular; GammaEtaBending; GammaEtaBendingRegular; GammaEtaBending2; GammaEtaBending2Regular; GammaYShearing; GammaYShearingRegular; GammaYShearing2; GammaYShearing2Regular; GammaXiBendingLinear; GammaEtaBendingLinear; GammaXiBendingLinear2; GammaEtaBendingLinear2; GammaXiBendingLinear3; GammaEtaBendingLinear3; alpha8; GammaXiBendingLinearRegular; GammaEtaBendingLinearRegular], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]) == subs([0; 0; 0; 0; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1;0;0;0;0;0;0;1/eta^2;0;0], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20]);
sol = solve([GammaXiBending; GammaXiBendingRegular; GammaXiBending2; GammaXiBending2Regular; GammaXShearing; GammaXShearingRegular; GammaXShearing2; GammaXShearing2Regular; GammaEtaBending; GammaEtaBendingRegular; GammaEtaBending2; GammaEtaBending2Regular; GammaYShearing; GammaYShearingRegular; GammaYShearing2; GammaYShearing2Regular; GammaXiBendingLinear; GammaEtaBendingLinear; GammaXiBendingLinear2; GammaEtaBendingLinear2; GammaXiBendingLinear3; GammaEtaBendingLinear3; alpha8; GammaXiBendingLinearRegular; GammaEtaBendingLinearRegular] == [0; 0; 0; 0; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1;0;0;0;0;0;0;1/eta^2;0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20]);
% sol = solve(subs([GammaXiBending; GammaEtaBending; GammaXiBendingRegular; GammaEtaBendingRegular; GammaXiBending2; GammaEtaBending2; GammaXiBending2Regular; GammaEtaBending2Regular; GammaXiShearing; GammaEtaShearing; GammaXiShearingLinear; GammaEtaShearingLinear], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]) == subs([0; 0; 0; 0; 0; 0; 0; 0; [J11, J21]*[1;0]; [J12, J22]*[1;0]; [J11, J21]*[y;x]; [J12, J22]*[y;x]], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20]);
% sol = solve(subs([GammaXiBending; GammaXiBendingRegular; GammaXiBending2; GammaXiBending2Regular; GammaXShearing; GammaXShearingRegular; GammaXShearing2; GammaXShearing2Regular; GammaXShearingLinearRegular; GammaEtaBending; GammaEtaBendingRegular; GammaEtaBending2; GammaEtaBending2Regular; GammaYShearing; GammaYShearingRegular; GammaYShearing2; GammaYShearing2Regular; GammaYShearingLinearRegular], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]) == subs([0; 0; 0; 0; 1; 1; 0; 0; subs(y, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]); 0; 0; 0; 0; 0; 0; 1; 1; subs(x, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1])], [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20]);
% solO = solve([GammaXiBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaEtaBending; GammaYBending2; GammaYShearing; GammaYShearing2; GammaXBending3; GammaYBending3] == [0; 0; 1; 0; 0; 0; 0; 1; 0; 0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
GammaXiCorrect = subs(subs(GammaXi, [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10]);
GammaEtaCorrect = subs(subs(GammaEta, [x1, x2, x3, x4, y1, y2, y3, y4], [1, 3,4,3/2,1,2,4,4]), [alpha11, alpha12, alpha13, alpha14, alpha15, alpha16, alpha17, alpha18, alpha19, alpha20], [sol.alpha11, sol.alpha12, sol.alpha13, sol.alpha14, sol.alpha15, sol.alpha16, sol.alpha17, sol.alpha18, sol.alpha19, sol.alpha20]);


GammaXCorrect = subs(GammaX, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]);
GammaYCorrect = subs(GammaY, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]);

ANSSolution_GammaXi = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
ANSSolution_GammaEta = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);
ANSSolution_GammaX = JInv(1, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];
ANSSolution_GammaY = JInv(2, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];

differenceGammaX = simplify(ANSSolution_GammaX - GammaXCorrect)
differenceGammaY = simplify(ANSSolution_GammaY - GammaYCorrect)

