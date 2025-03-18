%% Quadrilateral
syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12 alpha13 alpha14 alpha15 alpha16 alpha17 alpha18 alpha19 alpha20;
syms EsE;
syms C1 C2;

Es = [EsE, 0; 0, EsE];
Es = [1, 0; 0, 1];


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
detJ0 = det(J0);

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
% pointsInSkewCoordinates = [xi; eta] + (J0 \ c) * H1;
pointsInSkewCoordinates = [xi; eta] + [C1; C2] * H1;
nodalPointsInSkewCoordinates = cell(2, 4);
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates{1, ii} = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates{2, ii} = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
h = 1/4*[1, -1, 1, -1].';
C1test = 1/detJ0*(a2.'*yNodes.'*h.'*xNodes.'-a2.'*xNodes.'*h.'*yNodes.');
C2test = 1/detJ0*(-a1.'*yNodes.'*h.'*xNodes.'+a1.'*xNodes.'*h.'*yNodes.');

pointsInSkewCoordinates = [xiBar; etaBar];

% Metric shape functions
xPoints = [nodalPointsInSkewCoordinates{1, :}].';
yPoints = [nodalPointsInSkewCoordinates{2, :}].';
P = [ones(4, 1), xPoints, yPoints, xPoints .* yPoints];
p = [1, pointsInSkewCoordinates(1), pointsInSkewCoordinates(2), pointsInSkewCoordinates(1) .* pointsInSkewCoordinates(2)];
MInversion = p / P;
M_k_I =subs(MInversion, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);
M1 = M_k_I(1);
M2 = M_k_I(2);
M3 = M_k_I(3);
M4 = M_k_I(4);

M_k_I_Petrov = MInversion;
M1P = M_k_I_Petrov(1);
M2P = M_k_I_Petrov(2);
M3P = M_k_I_Petrov(3);
M4P = M_k_I_Petrov(4);

wPetrov = M1*w1+M2*w2+M3*w3+M4*w4;
thetaXPetrov = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
thetaYPetrov = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;

wPetrovBar = M1P*w1+M2P*w2+M3P*w3+M4P*w4;

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J011*thetaX+J021*thetaY;
GammaEtaO = diff(w, eta)+J012*thetaX+J022*thetaY;

GammaXiPetrovO = diff(wPetrov, xi)+J11*thetaXPetrov+J21*thetaYPetrov;
GammaEtaPetrovO = diff(wPetrov, eta)+J12*thetaXPetrov+J22*thetaYPetrov;

GammaXiBarO = subs(diff(wPetrovBar, xiBar)+J011*thetaXPetrov+J021*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);
GammaEtaBarO = subs(diff(wPetrovBar, etaBar)+J012*thetaXPetrov+J022*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);

ANSSolution_GammaXi = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
ANSSolution_GammaEta = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);

% test ANSSolution
wTest = xi^3*eta^4;
thetaXiTest = -diff(wTest, xi);
thetaEtaTest = -diff(wTest, eta);
wTest = 0;
thetaXiTest = eta;
thetaEtaTest = 0;
thetaXYTest = inv(J0.')*[thetaXiTest; thetaEtaTest];
thetaXTest = thetaXYTest(1);
thetaYTest = thetaXYTest(2);
wTestNodes = [subs(wTest, [xi, eta], [-1, -1]), subs(wTest, [xi, eta], [1, -1]), subs(wTest, [xi, eta], [1, 1]), subs(wTest, [xi, eta], [-1, 1])];
thetaXTestNodes = [subs(thetaXTest, [xi, eta], [-1, -1]), subs(thetaXTest, [xi, eta], [1, -1]), subs(thetaXTest, [xi, eta], [1, 1]), subs(thetaXTest, [xi, eta], [-1, 1])];
thetaYTestNodes = [subs(thetaYTest, [xi, eta], [-1, -1]), subs(thetaYTest, [xi, eta], [1, -1]), subs(thetaYTest, [xi, eta], [1, 1]), subs(thetaYTest, [xi, eta], [-1, 1])];
simplify(subs(ANSSolution_GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [wTestNodes, thetaXTestNodes, thetaYTestNodes]))
simplify(subs(ANSSolution_GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [wTestNodes, thetaXTestNodes, thetaYTestNodes]))



ANSSolution_GammaXiPetrov = 1/2*(1+eta)*subs(GammaXiPetrovO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiPetrovO, [xi, eta], [0, -1]);
ANSSolution_GammaEtaPetrov = 1/2*(1+xi)*subs(GammaEtaPetrovO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaPetrovO, [xi, eta], [-1, 0]);

[BsGammaXiOC, BsGammaXiOT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaOC, BsGammaEtaOT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsO = JInv*[BsGammaXiOC;BsGammaEtaOC];

[BsGammaXiANSC, BsGammaXiANST] = coeffs(ANSSolution_GammaXi, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaANSC, BsGammaEtaANST] = coeffs(ANSSolution_GammaEta, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsANS = JInv*[BsGammaXiANSC;BsGammaEtaANSC];

[BsGammaXiANSPetrovC, BsGammaXiANSPetrovT] = coeffs(ANSSolution_GammaXiPetrov, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaANSPetrovC, BsGammaEtaANSPetrovT] = coeffs(ANSSolution_GammaEtaPetrov, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsANSPetrov = JInv*[BsGammaXiANSPetrovC;BsGammaEtaANSPetrovC];

[BsGammaXiBarOC, BsGammaXiBarOT] = coeffs(GammaXiBarO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaBarOC, BsGammaEtaBarOT] = coeffs(GammaEtaBarO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsOPetrov = J0Inv*[BsGammaXiBarOC;BsGammaEtaBarOC];

tt1 = BsANS.'*Es*BsOPetrov;
tt2 = BsO.'*Es*BsANSPetrov;

% inttt2 = int(int(simplify(subs(subs(tt2, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1.2, 1, -1.1, -1.1, -1, 1.3, 1]), [C1, C2], [0, 0])), xi, [-1, 1]), [-1, 1]);

t1 = (-tt1+tt2);
t2 = t1*detJ;
% t3 = simplify(subs(t2, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1.2, 1, -1.1, -1.1, -1, 1.3, 1]));
t3 = simplify(subs(subs(t2, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1]), [C1, C2], [0, 0]));
t4 = int(t3(1:4, 1:4), xi, [-1, 1]);
t5 = int(t4, eta, [-1, 1]);
test1 = simplify(t5);


test2_1 = simplify(subs(subs(Es*BsOPetrov-J0*Es*[BsGammaXiBarOC;BsGammaEtaBarOC], [C1; C2], (J0 \ c)), [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1.2, 1, -1.1, -1.1, -1, 1.3, 1]));
test2_2 = simplify(subs(subs(Es*BsOPetrov-J0*Es*[BsGammaXiBarOC;BsGammaEtaBarOC], [C1; C2], (J0 \ c)), [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1]));

BsSearched = sym('BsSearched', [2, 12], 'real');
% Definition der Gleichung
equation = int(int((BsSearched.'*Es*BsO - BsO.'*Es*BsANS)*detJ, xi, -1, 1), eta, -1, 1) == 0;

eqr = int(BsSearched, xi, -1, 1) == 0;
% Lösung der Gleichung
solution = solve(eqr, BsSearched);
