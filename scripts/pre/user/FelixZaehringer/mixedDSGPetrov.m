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
syms A1 A2 A3 A4 B1 B2 B3 B4;

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
nodalPointsInSkewCoordinates = sym(zeros(2, 4));
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates(1, ii) = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates(2, ii) = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

pointsInSkewCoordinates = [xiBar; etaBar];

% Metric shape functions
xPoints = [nodalPointsInSkewCoordinates(1, :)].';
yPoints = [nodalPointsInSkewCoordinates(2, :)].';
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
thetaXPetrovBar = M1P*thetaX1+M2P*thetaX2+M3P*thetaX3+M4P*thetaX4;
thetaYPetrovBar = M1P*thetaY1+M2P*thetaY2+M3P*thetaY3+M4P*thetaY4;

GammaXiBar = subs(diff(wPetrovBar, xiBar)+J011*thetaXPetrov+J021*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);
GammaEtaBar = subs(diff(wPetrovBar, etaBar)+J012*thetaXPetrov+J022*thetaYPetrov, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);

[BsGammaXiBarC, BsGammaXiBarT] = coeffs(GammaXiBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaBarC, BsGammaEtaBarT] = coeffs(GammaEtaBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

Ls = [BsGammaXiBarC;BsGammaEtaBarC];

%% Integration in skew coordinates
% => Does not work: transverse shear locking
% intXiW = int(diff(wPetrovBar, xiBar), xiBar);
% intXiBar = int(J011*thetaXPetrovBar+J021*thetaYPetrovBar, xiBar);
% deltaGammaXi2 = subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 2)) - subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% deltaGammaXi3 = subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 3)) - subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% deltaGammaXi4 = subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 4)) - subs(wPetrovBar + intXiBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% 
% intEtaBar = int(J012*thetaXPetrovBar+J022*thetaYPetrovBar, etaBar);
% deltaGammaEta2 = subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 2)) - subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% deltaGammaEta3 = subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 3)) - subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% deltaGammaEta4 = subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 4)) - subs(wPetrovBar + intEtaBar, [xiBar; etaBar], nodalPointsInSkewCoordinates(:, 1));
% 
% % Modified transverse shear strains
% gammaXiBar = diff(M2P, xiBar)*deltaGammaXi2+diff(M3P, xiBar)*deltaGammaXi3+diff(M4P, xiBar)*deltaGammaXi4;
% gammaEtaBar = diff(M2P, etaBar)*deltaGammaEta2+diff(M3P, etaBar)*deltaGammaEta3+diff(M4P, etaBar)*deltaGammaEta4;
% 
% % Bs Matrices
% [BsGammaXiBarDSG, BsGammaXiBarDSGT] = coeffs(gammaXiBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
% [BsGammaEtaBarDSG, BsGammaEtaBarDSGT] = coeffs(gammaEtaBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
% 
% BsDSG = [BsGammaXiBarDSG; BsGammaEtaBarDSG];
% 
% BsDSGTest = subs(BsDSG, [C1, C2, J012, J021, htx, hty], [0, 0, 0, 0, 0, 0]);

%% Integration in natural coordinates
% Shear Gap Distributions
intXi = int(J11*thetaXPetrov+J21*thetaYPetrov, xi);
% deltaGammaXi2 = w2 - w1 + subs(intXi, [xi, eta], [1, -1]) - subs(intXi, [xi, eta], [-1, -1]);
% deltaGammaXi3 = w3 - w1 + subs(intXi, [xi, eta], [1, 1]) - subs(intXi, [xi, eta], [-1, -1]);
% deltaGammaXi4 = w4 - w1 + subs(intXi, [xi, eta], [-1, 1]) - subs(intXi, [xi, eta], [-1, -1]);

deltaGammaXi2 = w2 - w1 + subs(intXi, [xi, eta], [1, -1]) - subs(intXi, [xi, eta], [-1, -1]);
deltaGammaXi3 = w3 - w4 + subs(intXi, [xi, eta], [1, 1]) - subs(intXi, [xi, eta], [-1, 1]);
deltaGammaXi4 = 0;


% intXi = int(J11*thetaXPetrov+J21*thetaYPetrov, xi, [-1, 1]);
% deltaGammaXi2 = subs(intXi, [xi, eta], [1, -1]) - subs(intXi, [xi, eta], [-1, -1]);
% deltaGammaXi3 = subs(intXi, [xi, eta], [1, 1]) - subs(intXi, [xi, eta], [-1, -1]);
% deltaGammaXi4 = subs(intXi, [xi, eta], [-1, 1]) - subs(intXi, [xi, eta], [-1, -1]);

intEta = int(J12*thetaXPetrov+J22*thetaYPetrov, eta);
% deltaGammaEta2 = w2 - w1 + subs(intEta, [xi, eta], [1, -1]) - subs(intEta, [xi, eta], [-1, -1]);
% deltaGammaEta3 = w3 - w1 + subs(intEta, [xi, eta], [1, 1]) - subs(intEta, [xi, eta], [-1, -1]);
% deltaGammaEta4 = w4 - w1 + subs(intEta, [xi, eta], [-1, 1]) - subs(intEta, [xi, eta], [-1, -1]);

deltaGammaEta2 = 0;
deltaGammaEta3 = w3- w2 + subs(intEta, [xi, eta], [1, 1]) - subs(intEta, [xi, eta], [1, -1]);
deltaGammaEta4 = w4 - w1 + subs(intEta, [xi, eta], [-1, 1]) - subs(intEta, [xi, eta], [-1, -1]);


% Transformation
deltaGammaXiBar2 = subs(inv(J0)*J*[deltaGammaXi2;deltaGammaEta2], [xi, eta], [1, -1]);
deltaGammaXiBar3 = subs(inv(J0)*J*[deltaGammaXi3;deltaGammaEta3], [xi, eta], [1, 1]);
deltaGammaXiBar4 = subs(inv(J0)*J*[deltaGammaXi4;deltaGammaEta4], [xi, eta], [-1, 1]);

% Discrete Shear Gaps
% deltaGammaXiBar = M2P*deltaGammaXiBar2(1)+M3P*deltaGammaXiBar3(1)+M4P*deltaGammaXiBar4(1);
% deltaGammaEtaBar = M2P*deltaGammaXiBar2(2)+M3P*deltaGammaXiBar3(2)+M4P*deltaGammaXiBar4(2);

% Modified transverse shear strains
gammaXiBar = diff(M2P, xiBar)*deltaGammaXiBar2(1)+diff(M3P, xiBar)*deltaGammaXiBar3(1)+diff(M4P, xiBar)*deltaGammaXiBar4(1);
gammaEtaBar = diff(M2P, etaBar)*deltaGammaXiBar2(2)+diff(M3P, etaBar)*deltaGammaXiBar3(2)+diff(M4P, etaBar)*deltaGammaXiBar4(2);

% gammaXi = diff(M2, xi) * deltaGammaXi2 + diff(M3, xi) * deltaGammaXi3 + diff(M4, xi) * deltaGammaXi4;
% gammaEta = diff(M2, eta) * deltaGammaEta2 + diff(M3, eta) * deltaGammaEta3 + diff(M4, eta) * deltaGammaEta4;
% 
% gammaXiBar = diff(M1P, xiBar) * A1 + diff(M2P, xiBar) * A2 + diff(M3P, xiBar) * A3 + diff(M4P, xiBar) * A4;
% gammaEtaBar = diff(M1P, etaBar) * B1 + diff(M2P, etaBar) * B2 + diff(M3P, etaBar) * B3 + diff(M4P, etaBar) * B4;
% 
% gammaXiTrans = J.'*inv(J0.')*[gammaXiBar; gammaEtaBar];
% gammaXiTrans = subs(gammaXiTrans, [xiBar; etaBar], [xi; eta] + [C1; C2] * H1);
% 
% eq1 = subs(gammaXiTrans - [gammaXi; gammaEta], [xi, eta], [-1, -1]);
% eq2 = subs(gammaXiTrans - [gammaXi; gammaEta], [xi, eta], [1, -1]);
% eq3 = subs(gammaXiTrans - [gammaXi; gammaEta], [xi, eta], [1, 1]);
% eq4 = subs(gammaXiTrans - [gammaXi; gammaEta], [xi, eta], [-1, 1]);
% 
% sol = solve([eq1;eq2;eq3;eq4], [A1, A2, A3, A4, B1, B2, B3, B4]); % not solvable



% Bs Matrices
[BsGammaXiBarDSG, BsGammaXiBarDSGT] = coeffs(gammaXiBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[BsGammaEtaBarDSG, BsGammaEtaBarDSGT] = coeffs(gammaEtaBar, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BsDSG = [BsGammaXiBarDSG; BsGammaEtaBarDSG];

BsDSGTest = subs(BsDSG, [C1, C2, J012, J021, htx, hty], [0, 0, 0, 0, 0, 0]);

% test
int1 = int(diff(w, xi) + J11*thetaY+J21*thetaX, xi);
deltaGamma12 = subs(int1, [xi, eta], [1, -1]) - subs(int1, [xi, eta], [-1, -1]);
deltaGamma13 = subs(int1, [xi, eta], [1, 1]) - subs(int1, [xi, eta], [-1, -1]);
deltaGamma14 = subs(int1, [xi, eta], [-1, 1]) - subs(int1, [xi, eta], [-1, -1]);
GammaXiDSG = diff(N2, xi)*deltaGamma12+diff(N3, xi)*deltaGamma13+diff(N4, xi)*deltaGamma14;
[BsGammaXiDSG, BsGammaXiDSGT] = coeffs(GammaXiDSG, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
GammaXiDSGTest = subs(BsGammaXiDSG, [C1, C2, J012, J021, htx, hty], [0, 0, 0, 0, 0, 0])
