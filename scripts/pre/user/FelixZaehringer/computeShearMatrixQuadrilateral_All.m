% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms w1 w2 w3 w4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi omega MEI
syms xVar yVar
syms xiBar etaBar
syms etaBar1 etaBar2;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12 alpha13;

xNodes = [x1, x2, x3, x4];
yNodes = [y1, y2, y3, y4];

% Shape Functions
N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);

% Geometry
x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

x = subs(x, [x1 x2 x3 x4], xNodes);
y = subs(y, [y1 y2 y3 y4], yNodes);

% Jacobian
J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);
J = [J11, J12; J21, J22];

J0 = subs(J, {xi, eta}, {0, 0});
J011 = J0(1, 1);
J012 = J0(1, 2);
J021 = J0(2, 1);
J022 = J0(2, 2);
detJ0 = det(J0);

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

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;

ANSSolution_GammaXi = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
ANSSolution_GammaEta = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);

[ad, bd] = coeffs(subs(ANSSolution_GammaEta, [x1, x2, x3, x4, y1, y2, y3, y4], [1,8,9,2,1,3,16,14]), [xi, eta])

[dM1_xiBar, cXi1] = coeffs(subs(diff(M1, xiBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM2_xiBar, cXi2] = coeffs(subs(diff(M2, xiBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM3_xiBar, cXi3] = coeffs(subs(diff(M3, xiBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM4_xiBar, cXi4] = coeffs(subs(diff(M4, xiBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);

[dM1_etaBar, cEta1] = coeffs(subs(diff(M1, etaBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM2_etaBar, cEta2] = coeffs(subs(diff(M2, etaBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM3_etaBar, cEta3] = coeffs(subs(diff(M3, etaBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[dM4_etaBar, cEta4] = coeffs(subs(diff(M4, etaBar), [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);


% constant part Xi
alpha1 = dM1_xiBar(3);
alpha2 = dM2_xiBar(3);
alpha3 = dM3_xiBar(3);
alpha4 = dM4_xiBar(3);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaXi_Constant = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;J011;J021;J011;J021;0;0;0;J011*(y1+y2+y3+y4)/4-J021*(x1+x2+x3+x4)/4;0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% eta part Xi
alpha1 = dM1_xiBar(2);
alpha2 = dM2_xiBar(2);
alpha3 = dM3_xiBar(2);
alpha4 = dM4_xiBar(2);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaXi_Eta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J011*(y3/4 - y2/4 - y1/4 + y4/4)-J021*(x3/4 - x2/4 - x1/4 + x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% xi part Xi
alpha1 = 0;
alpha2 = 0;
alpha3 = 0;
alpha4 = 0;

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaXi_Xi = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J011*(y2/4 - y1/4 + y3/4 - y4/4)-J021*(x2/4 - x1/4 + x3/4 - x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% xi*eta part Xi
alpha1 = dM1_xiBar(1);
alpha2 = dM2_xiBar(1);
alpha3 = dM3_xiBar(1);
alpha4 = dM4_xiBar(1);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaXi_XiEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J011*(y1/4 - y2/4 + y3/4 - y4/4)-J021*(x1/4 - x2/4 + x3/4 - x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% constant part Eta
alpha1 = dM1_etaBar(3);
alpha2 = dM2_etaBar(3);
alpha3 = dM3_etaBar(3);
alpha4 = dM4_etaBar(3);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaEta_Constant = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;J012;J022;J012;J022;0;0;0;J012*(y1+y2+y3+y4)/4-J022*(x1+x2+x3+x4)/4;0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% eta part Eta
alpha1 = 0;
alpha2 = 0;
alpha3 = 0;
alpha4 = 0;

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaEta_Eta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J012*(y3/4 - y2/4 - y1/4 + y4/4)-J022*(x3/4 - x2/4 - x1/4 + x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% xi part Eta
alpha1 = dM1_etaBar(2);
alpha2 = dM2_etaBar(2);
alpha3 = dM3_etaBar(2);
alpha4 = dM4_etaBar(2);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaEta_Xi = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J012*(y2/4 - y1/4 + y3/4 - y4/4)-J022*(x2/4 - x1/4 + x3/4 - x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);

% xi*eta part Eta
alpha1 = dM1_etaBar(1);
alpha2 = dM2_etaBar(1);
alpha3 = dM3_etaBar(1);
alpha4 = dM4_etaBar(1);

eq1 = alpha1 +alpha2+alpha3+alpha4;
eq2 = alpha1*x1 + alpha2*x2 + alpha3*x3+ alpha4*x4;
eq3 = alpha1*y1 + alpha2*y2 + alpha3*y3+ alpha4*y4;
eq4 = alpha5+alpha6+alpha7+alpha8;
eq5 = alpha9+alpha10+alpha11+alpha12;
eq6 = alpha1*x1*y1 + alpha2*x2*y2 + alpha3*x3*y3+ alpha4*x4*y4 - alpha5*y1 - alpha6*y2-alpha7*y3-alpha8*y4-alpha9*x1-alpha10*x2-alpha11*x3-alpha12*x4;
eq7 = alpha1*x1^2 + alpha2*x2^2 + alpha3*x3^2+ alpha4*x4^2 - 2*alpha5*x1- 2*alpha6*x2-2*alpha7*x3- 2*alpha8*x4;
eq8 = alpha1*y1^2 + alpha2*y2^2 + alpha3*y3^2+ alpha4*y4^2 - 2*alpha9*y1- 2*alpha10*y2-2*alpha11*y3-2*alpha12*y4;
eq9 = alpha5*y1+alpha6*y2+alpha7*y3+alpha8*y4- alpha9*x1- alpha10*x2-alpha11*x3-alpha12*x4;
eq10 = alpha5-alpha6+alpha7-alpha8; % hourglass
eq11 = alpha9-alpha10+alpha11-alpha12; % hourglass

sol_GammaEta_XiEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11]==[0;0;0;0;0;0;0;0;J012*(y1/4 - y2/4 + y3/4 - y4/4)-J022*(x1/4 - x2/4 + x3/4 - x4/4);0;0], [alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);



% constant shearing part
eq1 = alpha1 + alpha2 +alpha3+alpha4;
eq2 = x1*alpha1 + x2*alpha2 +x3*alpha3+x4*alpha4;
eq3 = y1*alpha1 + y2*alpha2 +y3*alpha3+y4*alpha4;
eq4 = x1*y1*alpha1 + x2*y2*alpha2 +x3*y3*alpha3+x4*y4*alpha4;
% eq5 = 1/2*(x1^2*alpha1 + x2^2*alpha2 +x3^2*alpha3+x4^2*alpha4);

eq1 = 1/detJ0*(J022*(alpha1 + alpha2 +alpha3+alpha4)-J021*(alpha5 + alpha6 +alpha7+alpha8));
eq2 = 1/detJ0*(-J012*(alpha1 + alpha2 +alpha3+alpha4)+J011*(alpha5 + alpha6 +alpha7+alpha8));
eq3 = 1/detJ0*(J022*(x1*alpha1 + x2*alpha2 +x3*alpha3+x4*alpha4)-J021*(x1*alpha5 + x2*alpha6 +x3*alpha7+x4*alpha8));
eq4 = 1/detJ0*(-J012*(x1*alpha1 + x2*alpha2 +x3*alpha3+x4*alpha4)+J011*(x1*alpha5 + x2*alpha6 +x3*alpha7+x4*alpha8));
eq5 = 1/detJ0*(J022*(y1*alpha1 + y2*alpha2 +y3*alpha3+y4*alpha4)-J021*(y1*alpha5 + y2*alpha6 +y3*alpha7+y4*alpha8));
eq6 = 1/detJ0*(-J012*(y1*alpha1 + y2*alpha2 +y3*alpha3+y4*alpha4)+J011*(y1*alpha5 + y2*alpha6 +y3*alpha7+y4*alpha8));
eq7 = 1/detJ0*(J022*(x1*y1*alpha1 + x2*y2*alpha2 +x3*y3*alpha3+x4*y4*alpha4)-J021*(x1*y1*alpha5 + x2*y2*alpha6 +x3*y3*alpha7+x4*y4*alpha8));
eq8 = 1/detJ0*(-J012*(x1*y1*alpha1 + x2*y2*alpha2 +x3*y3*alpha3+x4*y4*alpha4)+J011*(x1*y1*alpha5 + x2*y2*alpha6 +x3*y3*alpha7+x4*y4*alpha8));

sol_w_ShearingConstant = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;1;0;0;1;(y1+y2+y3+y4)/4;(x1+x2+x3+x4)/4], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);

% sol_wxiBar_ShearingConstant = solve([eq1;eq2;eq3;eq4]==[0;J011;J021;J011*(y1+y2+y3+y4)/4+J021*(x1+x2+x3+x4)/4], [alpha1, alpha2, alpha3, alpha4]);
% sol_wetaBar_ShearingConstant = solve([eq1;eq2;eq3;eq4]==[0;J012;J022;J012*(y1+y2+y3+y4)/4+J022*(x1+x2+x3+x4)/4], [alpha1, alpha2, alpha3, alpha4]);

% eta shearing part

% eq1 = 1/detJ0*(J022*(alpha1 + alpha2 +alpha3+alpha4));
% eq2 = 1/detJ0*(-J012*(alpha1 + alpha2 +alpha3+alpha4));
% eq3 = 1/detJ0*(J022*(x1*alpha1 + x2*alpha2 +x3*alpha3+x4*alpha4));
% eq4 = 1/detJ0*(-J012*(x1*alpha1 + x2*alpha2 +x3*alpha3+x4*alpha4));
% eq5 = 1/detJ0*(J022*(y1*alpha1 + y2*alpha2 +y3*alpha3+y4*alpha4));
% eq6 = 1/detJ0*(-J012*(y1*alpha1 + y2*alpha2 +y3*alpha3+y4*alpha4));
% eq7 = 1/detJ0*(J022*(x1*y1*alpha1 + x2*y2*alpha2 +x3*y3*alpha3+x4*y4*alpha4));
% eq8 = 1/detJ0*(-J012*(x1*y1*alpha1 + x2*y2*alpha2 +x3*y3*alpha3+x4*y4*alpha4));

sol_w_ShearingEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;0;0;(y3/4 - y2/4 - y1/4 + y4/4);(x3/4 - x2/4 - x1/4 + x4/4)], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);

% sol_wxiBar_ShearingEta = solve([eq1;eq2;eq3;eq4]==[0;0;0;J011*(y3/4 - y2/4 - y1/4 + y4/4)+J021*(x3/4 - x2/4 - x1/4 + x4/4)], [alpha1, alpha2, alpha3, alpha4]);
% sol_wetaBar_ShearingEta = solve([eq1;eq2;eq3;eq4]==[0;0;0;J012*(y3/4 - y2/4 - y1/4 + y4/4)+J022*(x3/4 - x2/4 - x1/4 + x4/4)], [alpha1, alpha2, alpha3, alpha4]);

% xi shearing part
sol_w_ShearingXi = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;0;0;(y2/4 - y1/4 + y3/4 - y4/4);(x2/4 - x1/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
% sol_wxiBar_ShearingXi = solve([eq1;eq2;eq3;eq4]==[0;0;0;J011*(y2/4 - y1/4 + y3/4 - y4/4)+J021*(x2/4 - x1/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4]);
% sol_wetaBar_ShearingXi = solve([eq1;eq2;eq3;eq4]==[0;0;0;J012*(y2/4 - y1/4 + y3/4 - y4/4)+J022*(x2/4 - x1/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4]);

% xi*eta shearing part
sol_w_ShearingXiEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;0;0;(y1/4 - y2/4 + y3/4 - y4/4);(x1/4 - x2/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
% sol_wxiBar_ShearingXiEta = solve([eq1;eq2;eq3;eq4]==[0;0;0;J011*(y1/4 - y2/4 + y3/4 - y4/4)+J021*(x1/4 - x2/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4]);
% sol_wetaBar_ShearingXiEta = solve([eq1;eq2;eq3;eq4]==[0;0;0;J012*(y1/4 - y2/4 + y3/4 - y4/4)+J022*(x1/4 - x2/4 + x3/4 - x4/4)], [alpha1, alpha2, alpha3, alpha4]);

save solShearingBendingFactors20231012 sol_wxiBar_ShearingConstant sol_wxiBar_ShearingEta sol_wxiBar_ShearingXi sol_wxiBar_ShearingXiEta sol_wetaBar_ShearingConstant sol_wetaBar_ShearingEta sol_wetaBar_ShearingXi sol_wetaBar_ShearingXiEta;

% constant bending part Xi
eq1 = sol_wxiBar_ShearingConstant.alpha1*x1 + sol_wxiBar_ShearingConstant.alpha2*x2 +sol_wxiBar_ShearingConstant.alpha3*x3+sol_wxiBar_ShearingConstant.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wxiBar_ShearingConstant.alpha1*y1 + sol_wxiBar_ShearingConstant.alpha2*y2 +sol_wxiBar_ShearingConstant.alpha3*y3+sol_wxiBar_ShearingConstant.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wxiBar_ShearingConstant.alpha1*x1^2 + sol_wxiBar_ShearingConstant.alpha2*x2^2 +sol_wxiBar_ShearingConstant.alpha3*x3^2+sol_wxiBar_ShearingConstant.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wxiBar_ShearingConstant.alpha1*y1^2 + sol_wxiBar_ShearingConstant.alpha2*y2^2 +sol_wxiBar_ShearingConstant.alpha3*y3^2+sol_wxiBar_ShearingConstant.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wxiBar_ShearingConstant.alpha1*x1*y1 + sol_wxiBar_ShearingConstant.alpha2*x2*y2 +sol_wxiBar_ShearingConstant.alpha3*x3*y3+sol_wxiBar_ShearingConstant.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wxiBar_ShearingConstant.alpha1*x1^2*y1 + sol_wxiBar_ShearingConstant.alpha2*x2^2*y2 +sol_wxiBar_ShearingConstant.alpha3*x3^2*y3+sol_wxiBar_ShearingConstant.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wxiBar_ShearingConstant.alpha1*x1*y1^2 + sol_wxiBar_ShearingConstant.alpha2*x2*y2^2 +sol_wxiBar_ShearingConstant.alpha3*x3*y3^2+sol_wxiBar_ShearingConstant.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wxiBar_ShearingConstant.alpha1*x1^2*y1^2 + sol_wxiBar_ShearingConstant.alpha2*x2^2*y2^2 +sol_wxiBar_ShearingConstant.alpha3*x3^2*y3^2+sol_wxiBar_ShearingConstant.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;
% eq3 = alpha1*x1 + alpha2*x2 + alpha3*x3 + alpha4*x4;
% eq4 = alpha5*y1 + alpha6*y2 + alpha7*y3 + alpha8*y4;
eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
% eq7 = sol_wxiBar_ShearingConstant.alpha1*x1^2*y1 + sol_wxiBar_ShearingConstant.alpha2*x2^2*y2 +sol_wxiBar_ShearingConstant.alpha3*x3^2*y3+sol_wxiBar_ShearingConstant.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq8 = sol_wxiBar_ShearingConstant.alpha1*x1*y1^2 + sol_wxiBar_ShearingConstant.alpha2*x2*y2^2 +sol_wxiBar_ShearingConstant.alpha3*x3*y3^2+sol_wxiBar_ShearingConstant.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
eq7 = sol_wxiBar_ShearingConstant.alpha1*x1^2*y1 + sol_wxiBar_ShearingConstant.alpha2*x2^2*y2 +sol_wxiBar_ShearingConstant.alpha3*x3^2*y3+sol_wxiBar_ShearingConstant.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4;
eq8 = sol_wxiBar_ShearingConstant.alpha1*x1*y1^2 + sol_wxiBar_ShearingConstant.alpha2*x2*y2^2 +sol_wxiBar_ShearingConstant.alpha3*x3*y3^2+sol_wxiBar_ShearingConstant.alpha4*x4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;

% sol_GammaXi_BendingConstant = solve(subs([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;0;0;0;0], [x1, x2, x3, x4, y1, y2, y3, y4], [-1,1,1,-1,-1,-1,1,1]), [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);

sol_GammaXi_BendingConstant = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J011*(y1+y2+y3+y4)/4;J021*(x1+x2+x3+x4)/4;J021*((x1+x2+x3+x4)/4)^2;J011*((y1+y2+y3+y4)/4)^2], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaXi_BendingConstant');


% eta bending part Xi
eq1 = sol_wxiBar_ShearingEta.alpha1*x1 + sol_wxiBar_ShearingEta.alpha2*x2 +sol_wxiBar_ShearingEta.alpha3*x3+sol_wxiBar_ShearingEta.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wxiBar_ShearingEta.alpha1*y1 + sol_wxiBar_ShearingEta.alpha2*y2 +sol_wxiBar_ShearingEta.alpha3*y3+sol_wxiBar_ShearingEta.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wxiBar_ShearingEta.alpha1*x1^2 + sol_wxiBar_ShearingEta.alpha2*x2^2 +sol_wxiBar_ShearingEta.alpha3*x3^2+sol_wxiBar_ShearingEta.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wxiBar_ShearingEta.alpha1*y1^2 + sol_wxiBar_ShearingEta.alpha2*y2^2 +sol_wxiBar_ShearingEta.alpha3*y3^2+sol_wxiBar_ShearingEta.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wxiBar_ShearingEta.alpha1*x1*y1 + sol_wxiBar_ShearingEta.alpha2*x2*y2 +sol_wxiBar_ShearingEta.alpha3*x3*y3+sol_wxiBar_ShearingEta.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wxiBar_ShearingEta.alpha1*x1^2*y1 + sol_wxiBar_ShearingEta.alpha2*x2^2*y2 +sol_wxiBar_ShearingEta.alpha3*x3^2*y3+sol_wxiBar_ShearingEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wxiBar_ShearingEta.alpha1*x1*y1^2 + sol_wxiBar_ShearingEta.alpha2*x2*y2^2 +sol_wxiBar_ShearingEta.alpha3*x3*y3^2+sol_wxiBar_ShearingEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wxiBar_ShearingEta.alpha1*x1^2*y1^2 + sol_wxiBar_ShearingEta.alpha2*x2^2*y2^2 +sol_wxiBar_ShearingEta.alpha3*x3^2*y3^2+sol_wxiBar_ShearingEta.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;

eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wxiBar_ShearingEta.alpha1*x1^2*y1 + sol_wxiBar_ShearingEta.alpha2*x2^2*y2 +sol_wxiBar_ShearingEta.alpha3*x3^2*y3+sol_wxiBar_ShearingEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wxiBar_ShearingEta.alpha1*x1*y1^2 + sol_wxiBar_ShearingEta.alpha2*x2*y2^2 +sol_wxiBar_ShearingEta.alpha3*x3*y3^2+sol_wxiBar_ShearingEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;


sol_GammaXi_BendingEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J011*(y3/4 - y2/4 - y1/4 + y4/4);J021*(x3/4 - x2/4 - x1/4 + x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaXi_BendingEta');

% xi bending part Xi
eq1 = sol_wxiBar_ShearingXi.alpha1*x1 + sol_wxiBar_ShearingXi.alpha2*x2 +sol_wxiBar_ShearingXi.alpha3*x3+sol_wxiBar_ShearingXi.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wxiBar_ShearingXi.alpha1*y1 + sol_wxiBar_ShearingXi.alpha2*y2 +sol_wxiBar_ShearingXi.alpha3*y3+sol_wxiBar_ShearingXi.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wxiBar_ShearingXi.alpha1*x1^2 + sol_wxiBar_ShearingXi.alpha2*x2^2 +sol_wxiBar_ShearingXi.alpha3*x3^2+sol_wxiBar_ShearingXi.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wxiBar_ShearingXi.alpha1*y1^2 + sol_wxiBar_ShearingXi.alpha2*y2^2 +sol_wxiBar_ShearingXi.alpha3*y3^2+sol_wxiBar_ShearingXi.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wxiBar_ShearingXi.alpha1*x1*y1 + sol_wxiBar_ShearingXi.alpha2*x2*y2 +sol_wxiBar_ShearingXi.alpha3*x3*y3+sol_wxiBar_ShearingXi.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wxiBar_ShearingXi.alpha1*x1^2*y1 + sol_wxiBar_ShearingXi.alpha2*x2^2*y2 +sol_wxiBar_ShearingXi.alpha3*x3^2*y3+sol_wxiBar_ShearingXi.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wxiBar_ShearingXi.alpha1*x1*y1^2 + sol_wxiBar_ShearingXi.alpha2*x2*y2^2 +sol_wxiBar_ShearingXi.alpha3*x3*y3^2+sol_wxiBar_ShearingXi.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wxiBar_ShearingXi.alpha1*x1^2*y1^2 + sol_wxiBar_ShearingXi.alpha2*x2^2*y2^2 +sol_wxiBar_ShearingXi.alpha3*x3^2*y3^2+sol_wxiBar_ShearingXi.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;
eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wxiBar_ShearingXi.alpha1*x1^2*y1 + sol_wxiBar_ShearingXi.alpha2*x2^2*y2 +sol_wxiBar_ShearingXi.alpha3*x3^2*y3+sol_wxiBar_ShearingXi.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wxiBar_ShearingXi.alpha1*x1*y1^2 + sol_wxiBar_ShearingXi.alpha2*x2*y2^2 +sol_wxiBar_ShearingXi.alpha3*x3*y3^2+sol_wxiBar_ShearingXi.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;

sol_GammaXi_BendingXi = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J011*(y2/4 - y1/4 + y3/4 - y4/4);J021*(x2/4 - x1/4 + x3/4 - x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaXi_BendingXi');

% xi*eta bending part Xi
eq1 = sol_wxiBar_ShearingXiEta.alpha1*x1 + sol_wxiBar_ShearingXiEta.alpha2*x2 +sol_wxiBar_ShearingXiEta.alpha3*x3+sol_wxiBar_ShearingXiEta.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wxiBar_ShearingXiEta.alpha1*y1 + sol_wxiBar_ShearingXiEta.alpha2*y2 +sol_wxiBar_ShearingXiEta.alpha3*y3+sol_wxiBar_ShearingXiEta.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wxiBar_ShearingXiEta.alpha1*x1^2 + sol_wxiBar_ShearingXiEta.alpha2*x2^2 +sol_wxiBar_ShearingXiEta.alpha3*x3^2+sol_wxiBar_ShearingXiEta.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wxiBar_ShearingXiEta.alpha1*y1^2 + sol_wxiBar_ShearingXiEta.alpha2*y2^2 +sol_wxiBar_ShearingXiEta.alpha3*y3^2+sol_wxiBar_ShearingXiEta.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wxiBar_ShearingXiEta.alpha1*x1*y1 + sol_wxiBar_ShearingXiEta.alpha2*x2*y2 +sol_wxiBar_ShearingXiEta.alpha3*x3*y3+sol_wxiBar_ShearingXiEta.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wxiBar_ShearingXiEta.alpha1*x1^2*y1 + sol_wxiBar_ShearingXiEta.alpha2*x2^2*y2 +sol_wxiBar_ShearingXiEta.alpha3*x3^2*y3+sol_wxiBar_ShearingXiEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wxiBar_ShearingXiEta.alpha1*x1*y1^2 + sol_wxiBar_ShearingXiEta.alpha2*x2*y2^2 +sol_wxiBar_ShearingXiEta.alpha3*x3*y3^2+sol_wxiBar_ShearingXiEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wxiBar_ShearingXiEta.alpha1*x1^2*y1^2 + sol_wxiBar_ShearingXiEta.alpha2*x2^2*y2^2 +sol_wxiBar_ShearingXiEta.alpha3*x3^2*y3^2+sol_wxiBar_ShearingXiEta.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;
eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wxiBar_ShearingXiEta.alpha1*x1^2*y1 + sol_wxiBar_ShearingXiEta.alpha2*x2^2*y2 +sol_wxiBar_ShearingXiEta.alpha3*x3^2*y3+sol_wxiBar_ShearingXiEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wxiBar_ShearingXiEta.alpha1*x1*y1^2 + sol_wxiBar_ShearingXiEta.alpha2*x2*y2^2 +sol_wxiBar_ShearingXiEta.alpha3*x3*y3^2+sol_wxiBar_ShearingXiEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;


sol_GammaXi_BendingXiEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J011*(y1/4 - y2/4 + y3/4 - y4/4);J021*(x1/4 - x2/4 + x3/4 - x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaXi_BendingXiEta');

% constant bending part Eta
eq1 = sol_wetaBar_ShearingConstant.alpha1*x1 + sol_wetaBar_ShearingConstant.alpha2*x2 +sol_wetaBar_ShearingConstant.alpha3*x3+sol_wetaBar_ShearingConstant.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wetaBar_ShearingConstant.alpha1*y1 + sol_wetaBar_ShearingConstant.alpha2*y2 +sol_wetaBar_ShearingConstant.alpha3*y3+sol_wetaBar_ShearingConstant.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wetaBar_ShearingConstant.alpha1*x1^2 + sol_wetaBar_ShearingConstant.alpha2*x2^2 +sol_wetaBar_ShearingConstant.alpha3*x3^2+sol_wetaBar_ShearingConstant.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wetaBar_ShearingConstant.alpha1*y1^2 + sol_wetaBar_ShearingConstant.alpha2*y2^2 +sol_wetaBar_ShearingConstant.alpha3*y3^2+sol_wetaBar_ShearingConstant.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wetaBar_ShearingConstant.alpha1*x1*y1 + sol_wetaBar_ShearingConstant.alpha2*x2*y2 +sol_wetaBar_ShearingConstant.alpha3*x3*y3+sol_wetaBar_ShearingConstant.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wetaBar_ShearingConstant.alpha1*x1^2*y1 + sol_wetaBar_ShearingConstant.alpha2*x2^2*y2 +sol_wetaBar_ShearingConstant.alpha3*x3^2*y3+sol_wetaBar_ShearingConstant.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wetaBar_ShearingConstant.alpha1*x1*y1^2 + sol_wetaBar_ShearingConstant.alpha2*x2*y2^2 +sol_wetaBar_ShearingConstant.alpha3*x3*y3^2+sol_wetaBar_ShearingConstant.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wetaBar_ShearingConstant.alpha1*x1^2*y1^2 + sol_wetaBar_ShearingConstant.alpha2*x2^2*y2^2 +sol_wetaBar_ShearingConstant.alpha3*x3^2*y3^2+sol_wetaBar_ShearingConstant.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;

eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wetaBar_ShearingConstant.alpha1*x1^2*y1 + sol_wetaBar_ShearingConstant.alpha2*x2^2*y2 +sol_wetaBar_ShearingConstant.alpha3*x3^2*y3+sol_wetaBar_ShearingConstant.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wetaBar_ShearingConstant.alpha1*x1*y1^2 + sol_wetaBar_ShearingConstant.alpha2*x2*y2^2 +sol_wetaBar_ShearingConstant.alpha3*x3*y3^2+sol_wetaBar_ShearingConstant.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;


sol_GammaEta_BendingConstant = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J012*(y1+y2+y3+y4)/4;J022*(x1+x2+x3+x4)/4;0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaEta_BendingConstant');

% eta bending part Eta
eq1 = sol_wetaBar_ShearingEta.alpha1*x1 + sol_wetaBar_ShearingEta.alpha2*x2 +sol_wetaBar_ShearingEta.alpha3*x3+sol_wetaBar_ShearingEta.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wetaBar_ShearingEta.alpha1*y1 + sol_wetaBar_ShearingEta.alpha2*y2 +sol_wetaBar_ShearingEta.alpha3*y3+sol_wetaBar_ShearingEta.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wetaBar_ShearingEta.alpha1*x1^2 + sol_wetaBar_ShearingEta.alpha2*x2^2 +sol_wetaBar_ShearingEta.alpha3*x3^2+sol_wetaBar_ShearingEta.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wetaBar_ShearingEta.alpha1*y1^2 + sol_wetaBar_ShearingEta.alpha2*y2^2 +sol_wetaBar_ShearingEta.alpha3*y3^2+sol_wetaBar_ShearingEta.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wetaBar_ShearingEta.alpha1*x1*y1 + sol_wetaBar_ShearingEta.alpha2*x2*y2 +sol_wetaBar_ShearingEta.alpha3*x3*y3+sol_wetaBar_ShearingEta.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wetaBar_ShearingEta.alpha1*x1^2*y1 + sol_wetaBar_ShearingEta.alpha2*x2^2*y2 +sol_wetaBar_ShearingEta.alpha3*x3^2*y3+sol_wetaBar_ShearingEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wetaBar_ShearingEta.alpha1*x1*y1^2 + sol_wetaBar_ShearingEta.alpha2*x2*y2^2 +sol_wetaBar_ShearingEta.alpha3*x3*y3^2+sol_wetaBar_ShearingEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wetaBar_ShearingEta.alpha1*x1^2*y1^2 + sol_wetaBar_ShearingEta.alpha2*x2^2*y2^2 +sol_wetaBar_ShearingEta.alpha3*x3^2*y3^2+sol_wetaBar_ShearingEta.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;

eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wetaBar_ShearingEta.alpha1*x1^2*y1 + sol_wetaBar_ShearingEta.alpha2*x2^2*y2 +sol_wetaBar_ShearingEta.alpha3*x3^2*y3+sol_wetaBar_ShearingEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wetaBar_ShearingEta.alpha1*x1*y1^2 + sol_wetaBar_ShearingEta.alpha2*x2*y2^2 +sol_wetaBar_ShearingEta.alpha3*x3*y3^2+sol_wetaBar_ShearingEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;

sol_GammaEta_BendingEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J012*(y3/4 - y2/4 - y1/4 + y4/4);J022*(x3/4 - x2/4 - x1/4 + x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaEta_BendingEta');

% xi bending part Eta
eq1 = sol_wetaBar_ShearingXi.alpha1*x1 + sol_wetaBar_ShearingXi.alpha2*x2 +sol_wetaBar_ShearingXi.alpha3*x3+sol_wetaBar_ShearingXi.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wetaBar_ShearingXi.alpha1*y1 + sol_wetaBar_ShearingXi.alpha2*y2 +sol_wetaBar_ShearingXi.alpha3*y3+sol_wetaBar_ShearingXi.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wetaBar_ShearingXi.alpha1*x1^2 + sol_wetaBar_ShearingXi.alpha2*x2^2 +sol_wetaBar_ShearingXi.alpha3*x3^2+sol_wetaBar_ShearingXi.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wetaBar_ShearingXi.alpha1*y1^2 + sol_wetaBar_ShearingXi.alpha2*y2^2 +sol_wetaBar_ShearingXi.alpha3*y3^2+sol_wetaBar_ShearingXi.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wetaBar_ShearingXi.alpha1*x1*y1 + sol_wetaBar_ShearingXi.alpha2*x2*y2 +sol_wetaBar_ShearingXi.alpha3*x3*y3+sol_wetaBar_ShearingXi.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wetaBar_ShearingXi.alpha1*x1^2*y1 + sol_wetaBar_ShearingXi.alpha2*x2^2*y2 +sol_wetaBar_ShearingXi.alpha3*x3^2*y3+sol_wetaBar_ShearingXi.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wetaBar_ShearingXi.alpha1*x1*y1^2 + sol_wetaBar_ShearingXi.alpha2*x2*y2^2 +sol_wetaBar_ShearingXi.alpha3*x3*y3^2+sol_wetaBar_ShearingXi.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wetaBar_ShearingXi.alpha1*x1^2*y1^2 + sol_wetaBar_ShearingXi.alpha2*x2^2*y2^2 +sol_wetaBar_ShearingXi.alpha3*x3^2*y3^2+sol_wetaBar_ShearingXi.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;
eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wetaBar_ShearingXi.alpha1*x1^2*y1 + sol_wetaBar_ShearingXi.alpha2*x2^2*y2 +sol_wetaBar_ShearingXi.alpha3*x3^2*y3+sol_wetaBar_ShearingXi.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wetaBar_ShearingXi.alpha1*x1*y1^2 + sol_wetaBar_ShearingXi.alpha2*x2*y2^2 +sol_wetaBar_ShearingXi.alpha3*x3*y3^2+sol_wetaBar_ShearingXi.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;



sol_GammaEta_BendingXi = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J012*(y2/4 - y1/4 + y3/4 - y4/4);J022*(x2/4 - x1/4 + x3/4 - x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaEta_BendingXi');

% xi*eta bending part Eta
eq1 = sol_wetaBar_ShearingXiEta.alpha1*x1 + sol_wetaBar_ShearingXiEta.alpha2*x2 +sol_wetaBar_ShearingXiEta.alpha3*x3+sol_wetaBar_ShearingXiEta.alpha4*x4 - alpha1 - alpha2 - alpha3 - alpha4;
eq2 = sol_wetaBar_ShearingXiEta.alpha1*y1 + sol_wetaBar_ShearingXiEta.alpha2*y2 +sol_wetaBar_ShearingXiEta.alpha3*y3+sol_wetaBar_ShearingXiEta.alpha4*y4 - alpha5 - alpha6 - alpha7 - alpha8;
eq3 = sol_wetaBar_ShearingXiEta.alpha1*x1^2 + sol_wetaBar_ShearingXiEta.alpha2*x2^2 +sol_wetaBar_ShearingXiEta.alpha3*x3^2+sol_wetaBar_ShearingXiEta.alpha4*x4^2 - alpha1*2*x1 - alpha2*2*x2 - alpha3*2*x3 - alpha4*2*x4;
eq4 = sol_wetaBar_ShearingXiEta.alpha1*y1^2 + sol_wetaBar_ShearingXiEta.alpha2*y2^2 +sol_wetaBar_ShearingXiEta.alpha3*y3^2+sol_wetaBar_ShearingXiEta.alpha4*y4^2 - alpha5*2*y1 - alpha6*2*y2 - alpha7*2*y3 - alpha8*2*y4;
% eq5 = sol_wetaBar_ShearingXiEta.alpha1*x1*y1 + sol_wetaBar_ShearingXiEta.alpha2*x2*y2 +sol_wetaBar_ShearingXiEta.alpha3*x3*y3+sol_wetaBar_ShearingXiEta.alpha4*x4*y4 - alpha1*y1- alpha2*y2 - alpha3*y3 - alpha4*y4 - alpha5*x1 - alpha6*x2 - alpha7*x3 - alpha8*x4;
% eq6 = sol_wetaBar_ShearingXiEta.alpha1*x1^2*y1 + sol_wetaBar_ShearingXiEta.alpha2*x2^2*y2 +sol_wetaBar_ShearingXiEta.alpha3*x3^2*y3+sol_wetaBar_ShearingXiEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
% eq7 = sol_wetaBar_ShearingXiEta.alpha1*x1*y1^2 + sol_wetaBar_ShearingXiEta.alpha2*x2*y2^2 +sol_wetaBar_ShearingXiEta.alpha3*x3*y3^2+sol_wetaBar_ShearingXiEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;
% eq8 = sol_wetaBar_ShearingXiEta.alpha1*x1^2*y1^2 + sol_wetaBar_ShearingXiEta.alpha2*x2^2*y2^2 +sol_wetaBar_ShearingXiEta.alpha3*x3^2*y3^2+sol_wetaBar_ShearingXiEta.alpha4*x4^2*y4^2 - alpha1*2*x1*y1^2 - alpha2*2*x2*y2^2 - alpha3*2*x3*y3^2 - alpha4*2*x4*y4^2 - alpha5*2*x1^2*y1 - alpha6*2*x2^2*y2 - alpha7*2*x3^2*y3 - alpha8*2*x4^2*y4;
eq5 = alpha1*y1 + alpha2*y2 + alpha3*y3 + alpha4*y4;
eq6 = alpha5*x1 + alpha6*x2 + alpha7*x3 + alpha8*x4;
eq7 = sol_wetaBar_ShearingXiEta.alpha1*x1^2*y1 + sol_wetaBar_ShearingXiEta.alpha2*x2^2*y2 +sol_wetaBar_ShearingXiEta.alpha3*x3^2*y3+sol_wetaBar_ShearingXiEta.alpha4*x4^2*y4 - alpha1*2*x1*y1 - alpha2*2*x2*y2 - alpha3*2*x3*y3 - alpha4*2*x4*y4 - alpha5*x1^2 - alpha6*x2^2 - alpha7*x3^2 - alpha8*x4^2;
eq8 = sol_wetaBar_ShearingXiEta.alpha1*x1*y1^2 + sol_wetaBar_ShearingXiEta.alpha2*x2*y2^2 +sol_wetaBar_ShearingXiEta.alpha3*x3*y3^2+sol_wetaBar_ShearingXiEta.alpha4*x4*y4^2 - alpha1*y1^2 - alpha2*y2^2 - alpha3*y3^2 - alpha4*y4^2 - alpha5*2*x1*y1 - alpha6*2*x2*y2 - alpha7*2*x3*y3 - alpha8*2*x4*y4;


sol_GammaEta_BendingXiEta = solve([eq1;eq2;eq3;eq4;eq5;eq6;eq7;eq8]==[0;0;0;0;J012*(y1/4 - y2/4 + y3/4 - y4/4);J022*(x1/4 - x2/4 + x3/4 - x4/4);0;0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
disp('Done sol_GammaEta_BendingXiEta');

save solShearingBendingFactors sol_wxiBar_ShearingConstant sol_wxiBar_ShearingEta sol_wxiBar_ShearingXi sol_wxiBar_ShearingXiEta sol_wetaBar_ShearingConstant sol_wetaBar_ShearingEta sol_wetaBar_ShearingXi sol_wetaBar_ShearingXiEta sol_GammaXi_BendingConstant sol_GammaXi_BendingEta sol_GammaXi_BendingXi sol_GammaXi_BendingXiEta sol_GammaEta_BendingConstant sol_GammaEta_BendingEta sol_GammaEta_BendingXi sol_GammaEta_BendingXiEta;
