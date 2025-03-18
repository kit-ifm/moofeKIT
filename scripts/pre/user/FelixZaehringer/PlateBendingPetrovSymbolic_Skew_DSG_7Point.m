% Script to investigate transverse shear locking for Bubnov-Galerkin FE

% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms w1 w2 w3 w4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi omega MEI
syms xVar yVar
syms xiBar etaBar

% plateLength
plateLength = 50;

% Nodal coordinates
xNodes = [0, 3/5*plateLength, 2/5*plateLength, 0];
yNodes = [0, 0, 3/5*plateLength, plateLength];

% xNodes = [0, plateLength, 3/2*plateLength, 1/8*plateLength];
% yNodes = [0, 0, plateLength, plateLength];

%standard
% xNodes = [0, 1, 1, 0];
% yNodes = [0, 0, 1, 1];

% xNodes = [x1, x2, x3, x4];
% yNodes = [y1, y2, y3, y4];

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

% displacement
w = M1*w1 + M2*w2 + M3*w3 + M4*w4;

wAna = subs(-1/2*MEI*xVar^2+MEI*x1*xVar-phi*xVar+omega+phi*x1-1/2*MEI*x1^2, x1, xNodes(1));
wAna2 = -1/2*MEI*yVar^2+MEI*y1*yVar-phi*yVar+omega+phi*y1-1/2*MEI*y1^2;
wSubs = simplify(subs(w, [w1 w2 w3 w4], [subs(wAna, xVar, xNodes(1)) subs(wAna, xVar, xNodes(2)) subs(wAna, xVar, xNodes(3)) subs(wAna, xVar, xNodes(4))]));
% wSubs = simplify(subs(w, [w1 w2 w3 w4], [subs(wAna, xVar, x1) subs(wAna, xVar, x2) subs(wAna, xVar, x3) subs(wAna, xVar, x4)]));
wSubs2 = simplify(subs(w, [w1 w2 w3 w4], [subs(wAna2, yVar, yNodes(1)) subs(wAna2, yVar, yNodes(2)) subs(wAna2, yVar, yNodes(3)) subs(wAna2, yVar, yNodes(4))]));

% Rotations
thetaX = M1*thetaX1 + M2*thetaX2 + M3*thetaX3 + M4*thetaX4;
thetaY = M1*thetaY1 + M2*thetaY2 + M3*thetaY3 + M4*thetaY4;

thetaYAna = -(-MEI*xVar+MEI*x1-phi);
thetaXAna2 = -(-MEI*yVar+MEI*y1-phi);

thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [subs(thetaYAna, xVar, xNodes(1)) subs(thetaYAna, xVar, xNodes(2)) subs(thetaYAna, xVar, xNodes(3)) subs(thetaYAna, xVar, xNodes(4))]));

thetaXSubs2 = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [subs(thetaXAna2, yVar, yNodes(1)) subs(thetaXAna2, yVar, yNodes(2)) subs(thetaXAna2, yVar, yNodes(3)) subs(thetaXAna2, yVar, yNodes(4))]));
thetaYSubs2 = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [0 0 0 0]));

% Backtransformation
int1 = 1/2*subs(int(diff(subs(wSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), xi) + subs(J11*thetaYSubs+J21*thetaXSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), xi), x1, xNodes(1));
deltaGamma11 = simplify(subs(int1, [xi, eta], [-1, -1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma12 = simplify(subs(int1, [xi, eta], [1, -1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma13 = simplify(subs(int1, [xi, eta], [1, 1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma14 = simplify(subs(int1, [xi, eta], [-1, 1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma15 = simplify(subs(int1, [xi, eta], [0.5, -1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma16 = simplify(subs(int1, [xi, eta], [0, 1]) - subs(int1, [xi, eta], [-1, -1]));
deltaGamma17 = simplify(subs(int1, [xi, eta], [0, 0]) - subs(int1, [xi, eta], [-1, -1]));


int2 = 1/2*subs(int(diff(subs(wSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), eta) + subs(J22*thetaXSubs+J12*thetaYSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), eta), x1, xNodes(1));
deltaGamma22 = simplify(subs(int2, [xi, eta], [1, -1]) - subs(int2, [xi, eta], [-1, -1]));
deltaGamma23 = simplify(subs(int2, [xi, eta], [1, 1]) - subs(int2, [xi, eta], [-1, -1]));
deltaGamma24 = simplify(subs(int2, [xi, eta], [-1, 1]) - subs(int2, [xi, eta], [-1, -1]));
deltaGamma25 = simplify(subs(int2, [xi, eta], [-1, 0]) - subs(int2, [xi, eta], [-1, -1]));
deltaGamma26 = simplify(subs(int2, [xi, eta], [1, 0]) - subs(int2, [xi, eta], [-1, -1]));
deltaGamma27 = simplify(subs(int2, [xi, eta], [0, 0]) - subs(int2, [xi, eta], [-1, -1]));

GammaEtaBarDSG = diff(M2, etaBar)*deltaGamma22+diff(M3, etaBar)*deltaGamma23+diff(M4, etaBar)*deltaGamma24;
GammaEtaBarDSG = simplify(GammaEtaBarDSG);


GammaDSG = J0.'\[GammaXiBarDSG; GammaEtaBarDSG];
