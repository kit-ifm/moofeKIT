syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12;
syms Ex Ey;
syms J011 J022 J012 J021;
syms htx hty;
syms C1 C2;
syms D0 D1 D2 D3 D4 D5;

Es = [Ex, 0; 0, Ey];

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

% J0 = subs(J, {xi, eta}, {0, 0});
% J011 = J0(1, 1);
% J012 = J0(1, 2);
% J021 = J0(2, 1);
% J022 = J0(2, 2);
J0 = [J011, J012; J021, J022];
detJ0 = det(J0);
J0Inv = inv(J0.');
 
J = J0 + [htx*eta, htx*xi; hty*eta, hty*xi];

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J011*thetaX+J021*thetaY;
GammaEtaO = diff(w, eta)+J012*thetaX+J022*thetaY;

[BsGammaXiOC, BsGammaXiOT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaOC, BsGammaEtaOT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsO = [BsGammaXiOC;BsGammaEtaOC];



S = [dirac(xi)*dirac(1-eta), 0, dirac(xi)*dirac(1+eta), 0; 0, dirac(xi+1)*dirac(eta), 0, dirac(xi-1)*dirac(eta)];

EO = [1/2*(1+eta), 0, 1/2*(1-eta), 0; 0, 1/2*(1-xi), 0, 1/2*(1+xi)];
ETest = [1/2*(D0*eta+D1*xi+D2*xi*eta+D3*xi*eta^2+D4*eta^2+D5), 0, 1/2*(1-eta+D0*eta+D2*xi*eta), 0; 0, 1/2*(1-xi), 0, 1/2*(1+xi)];
E = [1/2*(1+eta-3/2*xi^2+3/2*xi^2*eta), 0, 1/2*(1-eta-3/2*xi^2-3/2*xi^2*eta), 0; 0, 1/2*(1-xi-3/2*eta^2-3/2*xi*eta^2), 0, 1/2*(1+xi-3/2*eta^2+3/2*xi*eta^2)];
EBar = [1/2*(1+eta+C2*xi*eta), 0, 1/2*(1-eta-C2*xi*eta), 0; 0, 1/2*(1-xi-C1*xi*eta), 0, 1/2*(1+xi+C1*xi*eta)];
EHat = [1/2*(1+eta+C2*xi*eta), 0, 1/2*(1-eta-C2*xi*eta), 0; 0, 1/2*(1-xi-C1*xi*eta), 0, 1/2*(1+xi+C1*xi*eta)];

% evaluation of (2) => *(1/2) due to evalution at borders of integral (property of dirac delta in Matlab)
test1 = int(int(E.'*S*detJ0, xi, [-1, 1]), eta, [-1, 1]);

test4 = int(int(S.'*J.'*inv(J0.')*A, xi, [-1, 1]), eta, [-1, 1]);


% evaluation of (1)
desiredSolA = simplify(1/detJ0*[J022, -J012]*Es*inv(J0.')*subs(EHat, [xi, eta], [0, 1]));
desiredSolB = simplify(1/detJ0*[-J021, J011]*Es*inv(J0.')*subs(EHat, [xi, eta], [-1, 0]));
desiredSolC = simplify(1/detJ0*[J022, -J012]*Es*inv(J0.')*subs(EHat, [xi, eta], [0, -1]));
desiredSolD = simplify(1/detJ0*[-J021, J011]*Es*inv(J0.')*subs(EHat, [xi, eta], [1, 0]));
test2A = int(int(E(:, 1).'*inv(J0)*Es*inv(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]);
test2B = int(int(E(:, 2).'*inv(J0)*Es*inv(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]);
test2C = int(int(E(:, 3).'*inv(J0)*Es*inv(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]);
test2D = int(int(E(:, 4).'*inv(J0)*Es*inv(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]);
% correctSolution = [(4*Ex*J022^2)/3 + (4*Ey*J012^2)/3, - Ey*J011*J012 - Ex*J021*J022, (2*Ex*J022^2)/3 + (2*Ey*J012^2)/3, - Ey*J011*J012 - Ex*J021*J022];
% sol = solve(simplify(test2([1,3])-correctSolution([1,3])), [D3, D4]);
% test2Eq = E.'*inv(J0)*Es*inv(J0.')*E*detJ0;

T = int(int(BsO.'*adjoint(J)*J0*S*2, xi, [-1, 1]), eta, [-1, 1]);
TO = int(int(BsO.'*detJ*S*2, xi, [-1, 1]), eta, [-1, 1]);

eq1 = 4/3*D0+4/9*D2*C2==0;
eq2 = 4/9*C2*D0+4/9*D2==-4/9*C2;

sol = solve([eq1;eq2], [D0, D2]);


% new try2

desiredFirstRow = simplify(int(int(EO(:, 1).'*adjoint(J)*Es*inv(J0.')*EO, xi, [-1, 1]), eta, [-1, 1]));
firstRow = simplify(int(int(ETest(:, 1).'*adjoint(J)*Es*inv(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]));
firstEntryDiff = simplify(firstRow(1)-desiredFirstRow(1));
secondEntryDiff = simplify(firstRow(2)-desiredFirstRow(2));
thirdEntryDiff = simplify(firstRow(3)-desiredFirstRow(3));
fourthEntryDiff = simplify(firstRow(4)-desiredFirstRow(4));
% simplify(subs(firstEntryDiff, [D0, D2], [sol.D0, sol.D2]))
% simplify(subs(thirdEntryDiff, [D0, D2], [sol.D0, sol.D2]))
% simplify(subs(fourthEntryDiff, [D0, D2], [sol.D0, sol.D2]))
firstEq = int(int((1+eta+C2*xi*eta)*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==0;
secondEq = int(int((1+eta+C2*xi*eta)*xi*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==-4/9*C2;
thirdEq = int(int((1-eta-C2*xi*eta)*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==0;
fourthEq = int(int((1-eta-C2*xi*eta)*xi*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==+4/9*C2;
fifthEq = int(int((1+xi+C1*xi*eta)*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==0;
sixthEq = int(int((1+xi+C1*xi*eta)*xi*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==+4/9*C1;
seventhEq = int(int((1-xi-C1*xi*eta)*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==0;
eigthEq = int(int((1-xi-C1*xi*eta)*xi*ETest(1, 1), xi, [-1, 1]), eta, [-1, 1])==-4/9*C1;

LGS = [firstEq; secondEq; thirdEq; fourthEq; fifthEq; sixthEq; seventhEq; eigthEq];

solve(LGS, [D0, D1, D2, D3, D4, D5])


A = sym('A', [2, 4]);
inv(A);


