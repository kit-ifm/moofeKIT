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

N = [1/2*(1+eta),0,1/2*(1-eta), 0; 0, 1/2*(1-xi),0,1/2*(1+xi)];
NBar = [dirac(xi)*dirac(1-eta),0,dirac(xi)*dirac(1+eta), 0; 0, dirac(eta)*dirac(1+xi),0, dirac(eta)*dirac(1-xi)];

E = int(int(N.'*NBar*detJ, xi, [-1, 1]), eta, [-1, 1]);
ET = int(int(NBar.'*N*detJ, xi, [-1, 1]), eta, [-1, 1]);

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;

[GammaXiOC, GammaXiOCT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[GammaEtaOC, GammaEtaOCT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

B = [GammaXiOC; GammaEtaOC];

C = int(int(B.'*NBar*detJ, xi, [-1, 1]), eta, [-1, 1]);

GammaTest = inv(E.')*C.'*[w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4].';

% Passt !
assert(simplify(GammaTest(1) - subs(GammaXiO, [xi, eta], [0, 1])) == 0);
assert(simplify(GammaTest(2) - subs(GammaEtaO, [xi, eta], [-1, 0])) == 0);
assert(simplify(GammaTest(3) - subs(GammaXiO, [xi, eta], [0, -1])) == 0);
assert(simplify(GammaTest(4) - subs(GammaEtaO, [xi, eta], [1, 0])) == 0);

syms Emod nu h;
Es = eye(2, 2);
Es = 5 / 6 * Emod / (2 * (1 + nu)) * h * Es;

H = int(int(N.'*Es*N*detJ, xi, [-1, 1]), eta, [-1, 1]);

NTest = [1,0,eta, 0; 0, 1,0,xi];
ETest = int(int(NTest.'*NBar*detJ, xi, [-1, 1]), eta, [-1, 1]);
HTest = int(int(NTest.'*Es*N*detJ, xi, [-1, 1]), eta, [-1, 1]);

lambdaTestTest = -inv(ETest)*HTest*inv(E.')*C.'*[w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4].';

lambdaTest = -inv(E)*H*inv(E.')*C.'*[w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4].';

assert(all(simplify(H*GammaTest+E*lambdaTest)==zeros(4,1)));

KsTest = C*inv(E)*H*inv(E.')*C.';

BBar = [1/2*((1+eta)*subs(GammaXiOC, [xi, eta], [0, 1]) + (1-eta)*subs(GammaXiOC, [xi, eta], [0, -1])); 1/2*((1+xi)*subs(GammaEtaOC, [xi, eta], [1, 0]) + (1-xi)*subs(GammaEtaOC, [xi, eta], [-1, 0]))];

subs(KsTest, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1])

% BBarTest = N*inv(E.')*C.';
BBarTest = C*inv(E)*N.';

BBarTestTest = C*inv(ETest)*NTest.';

testBBar = eval(subs(subs(BBarTest.', [x1, x2, x3, x4, y1, y2, y3, y4], [0, 25, 25, 0, 0, 0, 25, 25]), [xi, eta], [-0.5774, -0.5774]));
testBBarAns = eval(subs(subs(BBar, [x1, x2, x3, x4, y1, y2, y3, y4], [0, 25, 25, 0, 0, 0, 25, 25]), [xi, eta], [-0.5774, -0.5774]));



%% Test semi Petrov procedure
ETest = int(int(NBar.'*N, xi, [-1, 1]), eta, [-1, 1]);
DTest = int(int(NBar.'*B, xi, [-1, 1]), eta, [-1, 1]);
GammaTest = inv(ETest)*DTest;

ETilde = int(int(N.'*N, xi, [-1, 1]), eta, [-1, 1]);
CTest = int(int(B.'*N, xi, [-1, 1]), eta, [-1, 1]);

BBarTest = simplify(N*inv(ETilde.')*CTest.');
lambdaAnaly = 0;
