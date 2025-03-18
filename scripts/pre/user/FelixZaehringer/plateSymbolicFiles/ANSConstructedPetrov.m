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

wIso = N1*w1+N2*w2+N3*w3+N4*w4;
thetaXIso = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaYIso = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(wIso, xi)+J11*thetaXIso+J21*thetaYIso;
GammaEtaO = diff(wIso, eta)+J12*thetaXIso+J22*thetaYIso;

[GammaXiOC, GammaXiOCT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[GammaEtaOC, GammaEtaOCT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BIso = [GammaXiOC; GammaEtaOC];


N = [1/2*(1+eta),0,1/2*(1-eta), 0; 0, 1/2*(1-xi),0,1/2*(1+xi)];
% NBar = subs([dirac(xiBar)*dirac(1-etaBar),0,dirac(xiBar)*dirac(1+etaBar), 0; 0, dirac(etaBar)*dirac(1+xiBar),0, dirac(etaBar)*dirac(1-xiBar)], [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1);
NBar = [dirac(xi)*dirac(1-eta),0,dirac(xi)*dirac(1+eta), 0; 0, dirac(eta)*dirac(1+xi),0, dirac(eta)*dirac(1-xi)];


E = int(int(N.'*inv(J)*J0*NBar*detJ, xi, [-1, 1]), eta, [-1, 1]);
ETest = sym(zeros(4,4));
ETest(:, 1) = subs(N.'*inv(J)*J0(:, 1)*detJ, [xi, eta], [0, 1]);
% ET = int(int(NBar.'*N*detJ, xi, [-1, 1]), eta, [-1, 1]);

C = int(int((JInv * BIso).'*J0*NBar*detJ, xi, [-1, 1]), eta, [-1, 1]);
% CT = int(int(NBar.'*J0.'*J0Inv*BIso*detJ, xi, [-1, 1]), eta, [-1, 1]);

Bbar = N*inv(E.')*C.';

testBBar = eval(subs(subs(Bbar, [x1, x2, x3, x4, y1, y2, y3, y4], [0, 25, 25, 0, 0, 0, 25, 25]), [xi, eta], [-0.5774, -0.5774]));

[GammaXiBarOC, GammaXiBarOCT] = coeffs(GammaXiBarO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[GammaEtaBarOC, GammaEtaBarOCT] = coeffs(GammaEtaBarO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

B = [GammaXiBarOC; GammaEtaBarOC];



% D = int(int(NBar.'*JInv.'*J0Inv*subs(B, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1)*detJ, xi, [-1, 1]), eta, [-1, 1]);

GammaTest = inv(E.')*D*[w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4].';

% Passt !
assert(simplify(GammaTest(1) - subs(subs(GammaXiBarO, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), [xi, eta], [0, 1])) == 0);
assert(simplify(GammaTest(2) - subs(GammaEtaBarO, [xi, eta], [-1, 0])) == 0);
assert(simplify(GammaTest(3) - subs(GammaXiBarO, [xi, eta], [0, -1])) == 0);
assert(simplify(GammaTest(4) - subs(GammaEtaBarO, [xi, eta], [1, 0])) == 0);

syms Emod nu h;
Es = eye(2, 2);
Es = 5 / 6 * Emod / (2 * (1 + nu)) * h * Es;

H = int(int(N.'*Es*N*detJ, xi, [-1, 1]), eta, [-1, 1]);

lambdaTest = -inv(E)*H*inv(E.')*C.'*[w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4].';

assert(all(simplify(H*GammaTest+E*lambdaTest)==zeros(4,1)));

KsTest = C*inv(E)*H*inv(E.')*C.';

BBar = [1/2*((1+eta)*subs(GammaXiBarOC, [xi, eta], [0, 1]) + (1-eta)*subs(GammaXiBarOC, [xi, eta], [0, -1])); 1/2*((1+xi)*subs(GammaEtaBarOC, [xi, eta], [1, 0]) + (1-xi)*subs(GammaEtaBarOC, [xi, eta], [-1, 0]))];

subs(KsTest, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1])

BBarTest = N*inv(E.')*C.';
