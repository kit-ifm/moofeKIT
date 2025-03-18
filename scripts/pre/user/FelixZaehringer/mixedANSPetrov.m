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
syms D0 D1 D2 D3 D4 D5 D6 D7 D8;

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
J11 = J(1,1);
J12 = J(1,2);
J21 = J(2,1);
J22 = J(2,2);

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

% Bubnov / Lagrange
w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;

[BsGammaXiOC, BsGammaXiOT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaOC, BsGammaEtaOT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsO = [BsGammaXiOC;BsGammaEtaOC];
BsANS = [1/2*(1+eta)*subs(BsGammaXiOC, [xi, eta], [0, 1])+1/2*(1-eta)*subs(BsGammaXiOC, [xi, eta], [0, -1]);...
    1/2*(1+xi)*subs(BsGammaEtaOC, [xi, eta], [1, 0])+1/2*(1-xi)*subs(BsGammaEtaOC, [xi, eta], [-1, 0])];
BsANSPetrov = [1/2*(1+eta)*subs(BsGammaXiOC, [xi, eta], [0, 1])+1/2*(1-eta)*subs(BsGammaXiOC, [xi, eta], [0, -1]);...
    1/2*(1+xi)*subs(BsGammaEtaOC, [xi, eta], [1, 0])+1/2*(1-xi)*subs(BsGammaEtaOC, [xi, eta], [-1, 0])];

BsEps = collect(BsANS-BsO, [xi, eta]);
BsEpsPetrov = collect(BsANSPetrov-BsO, [xi, eta]);

% Petrov / Metric
wPetrov = M1*w1+M2*w2+M3*w3+M4*w4;
thetaXPetrov = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
thetaYPetrov = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;

wPetrovBar = M1P*w1+M2P*w2+M3P*w3+M4P*w4;

GammaXiBar = subs(diff(wPetrovBar, xiBar)+J011*thetaXPetrov+J021*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);
GammaEtaBar = subs(diff(wPetrovBar, etaBar)+J012*thetaXPetrov+J022*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);

[BsGammaXiBarC, BsGammaXiBarT] = coeffs(GammaXiBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaBarC, BsGammaEtaBarT] = coeffs(GammaEtaBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

Ls = [BsGammaXiBarC;BsGammaEtaBarC];

% Approximation Matrices
SANS = 2*[dirac(xi)*dirac(1-eta), 0, dirac(xi)*dirac(1+eta), 0; 0, dirac(xi+1)*dirac(eta), 0, dirac(xi-1)*dirac(eta)];
S = [dirac(1-eta), 0, dirac(1+eta), 0; 0, dirac(xi+1), 0, dirac(xi-1)];
P = S;

EO = [1/2*(1+eta), 0, 1/2*(1-eta), 0; 0, 1/2*(1-xi), 0, 1/2*(1+xi)];
EOPetrov = [1/2*(1+eta+C2*xi*eta), 0, 1/2*(1-eta-C2*xi*eta), 0; 0, 1/2*(1-xi-C1*xi*eta), 0, 1/2*(1+xi+C1*xi*eta)];
SStressPetrov = [1, 0, eta+C2*xi*eta, 0; 0, 1, 0, xi+C1*xi*eta];
SStressPetrov = SStressPetrov.*S;
SStress = [1, 0, eta, 0; 0, 1, 0, xi];
SStressEval = S + SANS;

sigmaBubnov = SANS.*EO;
sigmaPetrov = SANS.*EOPetrov;

% ETest = [1/2*(D0*eta+D1*xi+D2*xi*eta+D3*xi*eta^2+D4*eta^2+D5), 0, 1/2*(1-eta+D0*eta+D2*xi*eta), 0; 0, 1/2*(1-xi), 0, 1/2*(1+xi)];
ETest = [1/2*(1+eta+D1*eta), 0, 1/2*(1-eta), 0; 0, 1/2*(1-xi), 0, 1/2*(1+xi)];
E = [1/2*(1+eta-3/2*xi^2+3/2*xi^2*eta), 0, 1/2*(1-eta-3/2*xi^2-3/2*xi^2*eta), 0; 0, 1/2*(1-xi-3/2*eta^2-3/2*xi*eta^2), 0, 1/2*(1+xi-3/2*eta^2+3/2*xi*eta^2)];
EBar = [1/2*(1+eta+C2*xi*eta), 0, 1/2*(1-eta-C2*xi*eta), 0; 0, 1/2*(1-xi-C1*xi*eta), 0, 1/2*(1+xi+C1*xi*eta)];
EHat = [1/2*(1+eta+C2*xi*eta), 0, 1/2*(1-eta-C2*xi*eta), 0; 0, 1/2*(1-xi-C1*xi*eta), 0, 1/2*(1+xi+C1*xi*eta)];

% evaluation of (2) => *(1/2) due to evalution at borders of integral (property of dirac delta in Matlab)
test1 = int(int(E.'*S*detJ0, xi, [-1, 1]), eta, [-1, 1]);
int(int((BsANS-BsO).'*detJ*dirac(eta-1)*dirac(xi), xi, [-1, 1]), eta, [-1, 1])
int(int((BsANS-BsO).'*adjoint(J)*J0*dirac(eta-1)*dirac(xi), xi, [-1, 1]), eta, [-1, 1])

% first try BsANSNuevo
% BsANSNuevo = sym(zeros(2,12));
% BsANSNuevo1 = 1/2*(1+eta)*subs(BsO, eta, 1)+1/2*(1-eta)*subs(BsO, eta, -1);
% BsANSNuevo1_part2 = 1/2*(1-eta)*subs(BsO, eta, -1);
% BsANSNuevo2 = 1/2*(1+xi)*subs(BsO, xi, 1)+1/2*(1-xi)*subs(BsO, xi, -1);
% BsANSNuevo(1, :) = BsANSNuevo1(1, :);
% BsANSNuevo(2, :) = BsANSNuevo2(2, :);
% BsANSNuevoEps = (BsANSNuevo-BsO);
% int(int(BsANSNuevo1_part2.'*adjoint(J)*J0*dirac(eta-1), xi, [-1, 1]), eta, [-1, 1])
% int(int(BsO.'*adjoint(J)*J0*dirac(eta-1), xi, [-1, 1]), eta, [-1, 1])
% int(int((BsANSNuevo-BsO).'*adjoint(J)*J0*dirac(eta-1), xi, [-1, 1]), eta, [-1, 1])
% int(int((BsANSNuevo-BsO).'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1])
% int(int(BsANSNuevoEps(:,12).'*adjoint(J)*J0*sigmaPetrov(:,1), xi, [-1, 1]), eta, [-1, 1])
% collect(int((BsANSNuevo-BsO).'*adjoint(J)*J0*eta,  eta, [-1, 1]), [xi, eta])

% second try BsANSNuevo: works
BsOBar = J0.'*adjoint(J.')*BsO;
BsANSNuevoBar = sym(zeros(2,12));
BsANSNuevoBar1 = 1/2*(1+eta+C2*xi*eta)*subs(BsOBar, [xi, eta], [0, 1])+1/2*(1-eta-C2*xi*eta)*subs(BsOBar, [xi, eta], [0, -1]);
BsANSNuevoBar2 = 1/2*(1+xi+C1*xi*eta)*subs(BsOBar, [xi, eta], [1, 0])+1/2*(1-xi-C1*xi*eta)*subs(BsOBar, [xi, eta], [-1, 0]);
% BsANSNuevoBar1 = 1/2*(1+eta)*subs(BsOBar, eta, 1)+1/2*(1-eta)*subs(BsOBar, eta, -1);
% BsANSNuevoBar2 = 1/2*(1+xi)*subs(BsOBar, xi, 1)+1/2*(1-xi)*subs(BsOBar, xi, -1);
BsANSNuevoBar(1, :) = BsANSNuevoBar1(1, :);
BsANSNuevoBar(2, :) = BsANSNuevoBar2(2, :);
BsANSNuevo = J.'*inv(J0.')*BsANSNuevoBar;

t0a = int(int((BsANSNuevoBar).'*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
t1a = int(int((BsANSNuevo).'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
t1b = int(int((BsO).'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
simplify(t0a-t1b);

% third try BsANSNuevo: does not work
BsOTilde = adjoint(J.')*BsO;
BsANSNuevoTilde = sym(zeros(2,12));
BsANSNuevoTilde1 = 1/2*(1+eta)*subs(BsOTilde, eta, 1)+1/2*(1-eta)*subs(BsOTilde, eta, -1);
BsANSNuevoTilde2 = 1/2*(1+xi)*subs(BsOTilde, xi, 1)+1/2*(1-xi)*subs(BsOTilde, xi, -1);
BsANSNuevoTilde(1, :) = BsANSNuevoTilde1(1, :);
BsANSNuevoTilde(2, :) = BsANSNuevoTilde2(2, :);
t2a = int(int((BsANSNuevoTilde).'*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
t2b = int(int((BsO).'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
t2c = int(int((BsOTilde).'*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]); % equivalent to t2b
simplify(t2a-t2b);

% fourth try
BsOBarN = J0.'*inv(J.')*BsO;
BsANSNuevoBarN = sym(zeros(2,12));
BsANSNuevoBarN1 = 1/2*(1+eta)*subs(BsOBarN, eta, 1)+1/2*(1-eta)*subs(BsOBarN, eta, -1);
BsANSNuevoBarN2 = 1/2*(1+xi)*subs(BsOBarN, xi, 1)+1/2*(1-xi)*subs(BsOBarN, xi, -1);
BsANSNuevoBarN(1, :) = BsANSNuevoBarN1(1, :);
BsANSNuevoBarN(2, :) = BsANSNuevoBarN2(2, :);

t4a = int(int((BsANSNuevoBarN(:,1)).'*sigmaPetrov*detJ, xi, [-1, 1]), eta, [-1, 1]);
t4b = int(int((BsO).'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
simplify(t0a-t1b);

% test the JA Theorem
JA = subs(J, eta, 1);
tJA1 = J0.'*adjoint(JA.')*subs(BsO, eta, 1);

% last try: sigma = deltaSigma
BsOTilde = BsO;
BsANSNuevoTilde = sym(zeros(2,12));
BsANSNuevoTilde1 = 1/2*(1+eta)*subs(BsOTilde, [xi, eta], [0, 1])+1/2*(1-eta)*subs(BsOTilde, [xi, eta], [0, -1]);
BsANSNuevoTilde2 = 1/2*(1+xi)*subs(BsOTilde, [xi, eta], [1, 0])+1/2*(1-xi)*subs(BsOTilde, [xi, eta], [-1, 0]);
BsANSNuevoTilde(1, :) = BsANSNuevoTilde1(1, :);
BsANSNuevoTilde(2, :) = BsANSNuevoTilde2(2, :);
t2a = int(int((BsANSNuevoTilde).'*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
t2b = int(int((BsO).'*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
simplify(t2a-t2b);


% new test for orthogonality on edges
orthogonalityEdgesBubnov = int(int(BsEps.'*sigmaBubnov, xi, [-1, 1]), eta, [-1, 1]); % satisfied
orthogonalityEdgesPetrov = int(int(BsEpsPetrov.'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]); % not satisfied

orthogonalityEdgesPetrovRightSide = int(int(BsO.'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);
orthogonalityEdgesPetrovLeftSide = int(int(BsANSPetrov.'*adjoint(J)*J0*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1]);

% test4 = int(int(S.'*J.'*inv(J0.')*A, xi, [-1, 1]), eta, [-1, 1]);

BsANS = simplify(int(int(BsO.'*SANS, xi, [-1, 1]), eta, [-1, 1]));
BsSRI = subs(BsO, [xi, eta], [0, 0]);
BsAlt = simplify(int(int(Ls.'*P, xi, [-1, 1]), eta, [-1, 1]));

T = simplify(int(int(BsO.'*P, xi, [-1, 1]), eta, [-1, 1]));

I = simplify(int(int(EO.'*P, xi, [-1, 1]), eta, [-1, 1]));

G = int(int(S.'*J.'*adjoint(J0.')*EHat, xi, [-1, 1]), eta, [-1, 1]);

A = int(int(S.'*J.'*adjoint(J0.')*Ls, xi, [-1, 1]), eta, [-1, 1]);

A = sym('A', [2, 4]);

% test orthogonality
BsBar = EO*inv(I.')*T.';
ETest = [1/2*(1+eta+D3*xi^2), 0, 1/2*(1-eta+D1*xi^2), 0; 0, 1/2*(1-xi+D2*eta^2), 0, 1/2*(1+xi)];
BsBarTest = ETest*inv(I.')*T.';

% errorBubnov = int(int((BsO-BsBar).'*SStress, xi, [-1, 1]), eta, [-1, 1]);
errorPetrovONE = int((BsO-BsBarTest).'*adjoint(J)*J0*SStressPetrov,  eta, [-1, 1]);
errorPetrov = int(int((BsO-BsBarTest).'*adjoint(J)*J0*SStressPetrov, xi, [-1, 1]), eta, [-1, 1]); % *detJ0

errorBubnovEval = subs(errorBubnov, [x1, x2, x3, x4; y1, y2, y3, y4], [-1, 1.2, 1.4, -0.95; -1, -0.9, 1.2, 0.95]);
errorPetrovEval = subs(errorPetrov, [x1, x2, x3, x4; y1, y2, y3, y4], [-1, 1.2, 1.4, -0.95; -1, -0.9, 1.2, 0.95]);

int(int(BsO.'*P, xi, [-1, 1]), eta, [-1, 1])
int(int(BsBar.'*P, xi, [-1, 1]), eta, [-1, 1])

% patch test
int(int((BsO-BsBar).'*1, xi, [-1, 1]), eta, [-1, 1])



int(int((BsO-BsSRI).'*adjoint(J)*J0*dirac(xi)*dirac(eta), xi, [-1, 1]), eta, [-1, 1])


% Test with EAS
% EEAS = [xi+0*eta+0*xi*eta-3*xi*eta^2+0*eta*xi^2+0*xi^2+0*eta^2+0*xi^2*eta^2, 0, xi*eta-C2*eta*xi^2, 0; 0, eta-3*xi^2*eta, 0, xi*eta-C1*xi*eta^2];
% EEAS = [xi-D1*eta^2, 0, xi*eta-1/3*C2*eta, 0; 0, eta, 0, xi*eta-1/3*C1*xi];
% EEAS = [xi-3*xi*eta^2, 0, -1/3*C2*eta+C2*eta*xi^2, 0; 0, eta-3*xi^2*eta, 0, xi*eta-1/3*C1*xi];
% EEAS = [xi-xi*eta^2, 0, xi*eta-1/3*C2*eta, 0; 0, eta, 0, xi*eta-1/3*C1*xi];
EEAS = [xi-1/2*C2*xi^2, 0, xi*eta-1/6*C2*eta, 0, xi*eta^2-xi, 0, eta^2-1, 0; 0, eta-1/2*C1*eta^2, 0, xi*eta-1/6*C1*xi, 0, xi^2*eta-eta, 0, xi^2-1];

% EEAS = [xi-3*xi*eta^2, 0, xi*eta+0*eta+D2*xi-3*D2*xi*eta^2+D4*eta*xi^2+D5*xi^2+D6*eta^2+D7*xi^2*eta^2, 0; 0, eta-3*xi^2*eta, 0, xi*eta-C1*xi*eta^2];
% EEAS = [xi+D1*eta+D2*xi*eta+D3*xi*eta^2+D4*eta*xi^2+D5*xi^2+D6*eta^2+D7*xi^2*eta^2, 0, xi*eta-C2*eta*xi^2, 0; 0, eta, 0, xi*eta-C1*xi*eta^2];


% int(int(EEAS.'*SStressPetrov, xi, [-1, 1]), eta, [-1, 1])
% 
% int(int(EEAS.'*inv(J0)*J*SANS, xi, [-1, 1]), eta, [-1, 1])

int(int(EEAS.'*sigmaPetrov, xi, [-1, 1]), eta, [-1, 1])

int(int(EEAS.'*inv(J0)*J*S, xi, [-1, 1]), eta, [-1, 1])
