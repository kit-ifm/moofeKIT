%% Quadrilateral
syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12;

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
nodalPointsInSkewCoordinates = sym(zeros(2, 4));
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates(1, ii) = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates(2, ii) = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

pointsInSkewCoordinates = [xiBar; etaBar];

% Metric shape functions
continuumObject.elementGeometryType = 'quadrilateral';
[M_k_I, ~] = computeMetricShapeFunctions(continuumObject, 2, nodalPointsInSkewCoordinates, pointsInSkewCoordinates);
M1 = M_k_I(1);
M2 = M_k_I(2);
M3 = M_k_I(3);
M4 = M_k_I(4);

w = M1*w1+M2*w2+M3*w3+M4*w4;
thetaX = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
thetaY = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;

GammaXiBarO = diff(w, xiBar)+J011*thetaX+J021*thetaY;
GammaEtaBarO = diff(w, etaBar)+J012*thetaX+J022*thetaY;

% Bbar
GammaXiBarEdt = subs(GammaXiBarO, xiBar, 0);
GammaEtaBarEdt = subs(GammaEtaBarO, etaBar, 0);

GammaXiBarDiff = subs(GammaXiBarO-GammaXiBarEdt, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1);

[GammaXiBarOC, GammaXiBarOCT] = coeffs(subs(GammaXiBarO, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta]);
[GammaXiBarDiffC, GammaXiBarDiffCT] = coeffs(GammaXiBarDiff, [xi, eta]);

% taken from ...
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

eqGammaXiBar = alpha1*w1 +alpha2*w2+alpha3*w3+alpha4*w4+alpha5*thetaX1+alpha6*thetaX2+alpha7*thetaX3+alpha8*thetaX4+alpha9*thetaY1+alpha10*thetaY2+alpha11*thetaY3+alpha12*thetaY4;

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

[GammaXiBarDiffC1, GammaXiBarDiff1CT] = coeffs(GammaXiBarDiff(1), [w1 w2 w3 w4 thetaX1 thetaX2 thetaX3 thetaX4 thetaY1 thetaY2 thetaY3 thetaY4]);
[GammaXiBarOC7, GammaXiBarOC7CT] = coeffs(GammaXiBarOC(7), [w1 w2 w3 w4 thetaX1 thetaX2 thetaX3 thetaX4 thetaY1 thetaY2 thetaY3 thetaY4]);

% int(int(GammaXiBarOC7(5)-alpha5, xi, [-1, 1]), eta, [-1, 1]);

eq9 = int(int((GammaXiBarDiffC1(1)-alpha5)*xi^2*eta^2+GammaXiBarOC7(5)-sol_GammaXi_Constant.alpha5, xi, [-1, 1]), eta, [-1, 1]);
alpha5Val = solve(eq9==0, alpha5); 
subs(alpha8Val, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,2,1.5,-1,-1.3,-1,1,1])

eq9 = int(int((GammaXiBarDiffC1(2)-alpha6)*xi^2*eta^2+GammaXiBarOC7(6)-sol_GammaXi_Constant.alpha6, xi, [-1, 1]), eta, [-1, 1]);
alpha6Val = solve(eq9==0, alpha6);

eq9 = int(int((GammaXiBarDiffC1(3)-alpha7)*xi^2*eta^2+GammaXiBarOC7(7)-sol_GammaXi_Constant.alpha7, xi, [-1, 1]), eta, [-1, 1]);
alpha7Val = solve(eq9==0, alpha7);

eq9 = int(int((GammaXiBarDiffC1(4)-alpha8)*xi^2*eta^2+GammaXiBarOC7(8)-sol_GammaXi_Constant.alpha8, xi, [-1, 1]), eta, [-1, 1]);
alpha8Val = solve(eq9==0, alpha8);

subs(alpha5Val+ alpha6Val+alpha7Val+alpha8Val, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,2,1.5,-0.8,-0.8,-1,1,1])

% Conclusion:
% Bbar müsste zusätzlich einen xi^2*eta^2-Term besitzen, damit Eigenschaft
% int(B-BBar)=0 eingehalten werden kann. Es ist allerdings mehr als
% fraglich, ob dqas funktionieren kann... => Nope, geht nicht...


% Nächster Schritt:
% 

