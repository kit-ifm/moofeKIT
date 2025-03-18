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
xNodes = [10, 3/5*plateLength, 2/5*plateLength, 0];
yNodes = [0, 0, 3/5*plateLength, plateLength];

% xNodes = [0, plateLength, 3/2*plateLength, 1/8*plateLength];
% yNodes = [0, 0, plateLength, plateLength];

%standard
% xNodes = [0, 1, 1, 0];
% yNodes = [0, 0, 1, 1];

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

JA = subs(J, {xi, eta}, {0, 1});
JA11 = JA(1, 1);
JA12 = JA(1, 2);
JA21 = JA(2, 1);
JA22 = JA(2, 2);
detJA = det(JA);

JB = subs(J, {xi, eta}, {-1, 0});
JB11 = JB(1, 1);
JB12 = JB(1, 2);
JB21 = JB(2, 1);
JB22 = JB(2, 2);
detJB = det(JB);

JC = subs(J, {xi, eta}, {0, -1});
JC11 = JC(1, 1);
JC12 = JC(1, 2);
JC21 = JC(2, 1);
JC22 = JC(2, 2);
detJC = det(JC);

JD = subs(J, {xi, eta}, {1, 0});
JD11 = JD(1, 1);
JD12 = JD(1, 2);
JD21 = JD(2, 1);
JD22 = JD(2, 2);
detJD = det(JD);

test1 = 1/detJB*(J012*JB22-J022*JB12)+1/detJD*(J012*JD22-J022*JD12);
test2 = 1/detJA*(J021*JA11-J011*JA21)+1/detJC*(J021*JC11-J011*JC21);

% skew coordinates
numberOfNodes = 4;
numberOfNodes_9 = 9;
nodesXi = elementNodesInLocalCoordinates(2, 'quadrilateral', 4);
nodesXi_9 = elementNodesInLocalCoordinates(2, 'quadrilateral', 9);

NSymbolic_1 = 1 / 4 * (xi.^2 - xi) .* (eta.^2 - eta);
NSymbolic_2 = 1 / 4 * (xi.^2 + xi) .* (eta.^2 - eta);
NSymbolic_3 = 1 / 4 * (xi.^2 + xi) .* (eta.^2 + eta);
NSymbolic_4 = 1 / 4 * (xi.^2 - xi) .* (eta.^2 + eta);
NSymbolic_5 = 1 / 2 * (1 - xi.^2) .* eta .* (eta - 1);
NSymbolic_6 = 1 / 2 * xi .* (xi + 1) .* (1 - eta.^2);
NSymbolic_7 = 1 / 2 * (1 - xi.^2) .* eta .* (eta + 1);
NSymbolic_8 = 1 / 2 * xi .* (xi - 1) .* (1 - eta.^2);
NSymbolic_9 = (1 - xi.^2) .* (1 - eta.^2);

NSymbolic = [NSymbolic_1; NSymbolic_2; NSymbolic_3; NSymbolic_4; NSymbolic_5; NSymbolic_6; NSymbolic_7; NSymbolic_8; NSymbolic_9];

NSymbolic_0 = subs(NSymbolic, [xi, eta], [0,0]);

xNodes_9 = [xNodes, 1/2*(xNodes(1)+xNodes(2)), 1/2*(xNodes(2)+xNodes(3)), 1/2*(xNodes(3)+xNodes(4)), 1/2*(xNodes(1)+xNodes(4)), 1/2*(1/2*(xNodes(1)+xNodes(4))+1/2*(xNodes(2)+xNodes(3)))];
yNodes_9 = [yNodes, 1/2*(yNodes(1)+yNodes(2)), 1/2*(yNodes(2)+yNodes(3)), 1/2*(yNodes(3)+yNodes(4)), 1/2*(yNodes(1)+yNodes(4)), 1/2*(1/2*(yNodes(1)+yNodes(4))+1/2*(yNodes(2)+yNodes(3)))];
x_9 = xNodes_9*NSymbolic;
y_9 = yNodes_9*NSymbolic;

J11_9 = diff(x_9, xi);
J12_9 = diff(x_9, eta);
J21_9 = diff(y_9, xi);
J22_9 = diff(y_9, eta);
J_9 = [J11_9, J12_9; J21_9, J22_9];

J0_9 = subs(J_9, {xi, eta}, {0, 0});

evaluationPointsInSkewCoordinates = J0_9 \ (nodesXi_9 * (NSymbolic - NSymbolic_0));

nodalPointsInSkewCoordinates_9 = cell(2, 9);
for ii = 1:numberOfNodes_9
    nodalPointsInSkewCoordinates_9{1, ii} = simplify(subs(evaluationPointsInSkewCoordinates(1), {xi; eta}, nodesXi_9(:, ii)));
    nodalPointsInSkewCoordinates_9{2, ii} = simplify(subs(evaluationPointsInSkewCoordinates(2), {xi; eta}, nodesXi_9(:, ii)));
end

evaluationPointsInSkewCoordinates = [xiBar; etaBar];

xPoints_9 = [nodalPointsInSkewCoordinates_9{1, :}].';
yPoints_9 = [nodalPointsInSkewCoordinates_9{2, :}].';
P_9 = [ones(9, 1), xPoints_9, yPoints_9, xPoints_9.^2, xPoints_9 .* yPoints_9, yPoints_9.^2, xPoints_9.^2 .* yPoints_9, xPoints_9 .* yPoints_9.^2, xPoints_9.^2 .* yPoints_9.^2];
p_9 = [1, evaluationPointsInSkewCoordinates(1), evaluationPointsInSkewCoordinates(2), evaluationPointsInSkewCoordinates(1)^2, evaluationPointsInSkewCoordinates(1) .* evaluationPointsInSkewCoordinates(2), evaluationPointsInSkewCoordinates(2)^2, evaluationPointsInSkewCoordinates(1)^2*evaluationPointsInSkewCoordinates(2), evaluationPointsInSkewCoordinates(1)*evaluationPointsInSkewCoordinates(2)^2, evaluationPointsInSkewCoordinates(1)^2*evaluationPointsInSkewCoordinates(2)^2];
M_k_I_9 = p_9 / P_9;

% check differences
% for ii=1:9
%     assert(subs(diff(M_k_I_9(ii), xiBar), [xiBar, etaBar], [0, 0]) == 0, ['error in', num2str(ii)]);
% end

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

w1Subs = simplify(subs(subs(wSubs, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]), [x1 x2 x3 x4], xNodes));
w2Subs = simplify(subs(subs(wSubs, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 2}]), [x1 x2 x3 x4], xNodes));
w3Subs = simplify(subs(subs(wSubs, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 3}]), [x1 x2 x3 x4], xNodes));
w4Subs = simplify(subs(subs(wSubs, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 4}]), [x1 x2 x3 x4], xNodes));

% Rotations
thetaX = M1*thetaX1 + M2*thetaX2 + M3*thetaX3 + M4*thetaX4;
thetaY = M1*thetaY1 + M2*thetaY2 + M3*thetaY3 + M4*thetaY4;

thetaYAna = -(-MEI*xVar+MEI*x1-phi);
thetaXAna2 = -(-MEI*yVar+MEI*y1-phi);

thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [subs(thetaYAna, xVar, xNodes(1)) subs(thetaYAna, xVar, xNodes(2)) subs(thetaYAna, xVar, xNodes(3)) subs(thetaYAna, xVar, xNodes(4))]));
thetaYSubs = simplify(subs(thetaYSubs, [x1 x2 x3 x4], xNodes));

thetaXSubs2 = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [subs(thetaXAna2, yVar, yNodes(1)) subs(thetaXAna2, yVar, yNodes(2)) subs(thetaXAna2, yVar, yNodes(3)) subs(thetaXAna2, yVar, yNodes(4))]));
thetaYSubs2 = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [0 0 0 0]));


% Transverse shear strain
% syms J011 J021 J012 J022
GammaXiBar = diff(wSubs, xiBar) + J011*thetaYSubs+J021*thetaXSubs;
GammaEtaBar = diff(wSubs, etaBar) + J022*thetaXSubs+J012*thetaYSubs;

% GammaXiBar = subs(GammaXiBar, x1, xNodes(1));

GammaXiBar2 = diff(wSubs2, xiBar) + J011*thetaYSubs2+J021*thetaXSubs2;
GammaEtaBar2 = diff(wSubs2, etaBar) + J022*thetaXSubs2+J012*thetaYSubs2;1

% GammaXiBar = simplify(collect(simplify(collect(collect(collect(GammaXiBar, omega), MEI), phi)), [xiBar, etaBar]));
GammaEtaBar2 = simplify(collect(simplify(collect(collect(collect(GammaEtaBar2, omega), MEI), phi)), [xiBar, etaBar]));

GammaXiBarANS = (1-etaBar)*subs(GammaXiBar, [xiBar, etaBar], [1/2*(nodalPointsInSkewCoordinates{1, 1} + nodalPointsInSkewCoordinates{1, 2}), 1/2*(nodalPointsInSkewCoordinates{2, 1} + nodalPointsInSkewCoordinates{2, 2})])+(1+etaBar)*subs(GammaXiBar, [xiBar, etaBar], [1/2*(nodalPointsInSkewCoordinates{1, 4} + nodalPointsInSkewCoordinates{1, 3}), 1/2*(nodalPointsInSkewCoordinates{2, 4} + nodalPointsInSkewCoordinates{2, 3})]);
pointA=solve(subs(GammaXiBar, etaBar, 1), xiBar);
pointB=solve(subs(GammaEtaBar2, xiBar, -1), etaBar);
pointC = solve(subs(GammaXiBar, etaBar, -1), xiBar);
pointD=solve(subs(GammaEtaBar2, xiBar, 1), etaBar);

% Backtransformation 
% Gamma = J0.'\[GammaXiBar; GammaEtaBar];
int1 = 1/2*subs(int(diff(wSubs, xiBar) + J011*thetaYSubs+J021*thetaXSubs, xiBar), x1, xNodes(1));
deltaGamma11 = simplify(subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]) - subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
deltaGamma12 = simplify(subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 2}]) - subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
deltaGamma13 = simplify(subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 3}]) - subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
deltaGamma14 = simplify(subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 4}]) - subs(int1, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
% 
int1Alt = int(diff(subs(wSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), xi) + subs(J11*thetaYSubs+J21*thetaXSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), xi);
deltaGamma12Alt = subs(int1Alt, [xi, eta], [1, -1]) - subs(int1Alt, [xi, eta], [-1, -1]);
deltaGamma13Alt = subs(int1Alt, [xi, eta], [1, 1]) - subs(int1Alt, [xi, eta], [-1, -1]);
deltaGamma14Alt = subs(int1Alt, [xi, eta], [-1, 1]) - subs(int1Alt, [xi, eta], [-1, -1]);
GammaXiBarDSG = diff(M2, xiBar)*deltaGamma12+diff(M3, xiBar)*deltaGamma13+diff(M4, xiBar)*deltaGamma14;
GammaXiBarDSG = simplify(GammaXiBarDSG);

int2 = 1/2*subs(int(diff(wSubs, etaBar) + J022*thetaXSubs+J012*thetaYSubs, etaBar), x1, xNodes(1));
deltaGamma22 = simplify(subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 2}]) - subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
deltaGamma23 = simplify(subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 3}]) - subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));
deltaGamma24 = simplify(subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 4}]) - subs(int2, [xiBar, etaBar], [nodalPointsInSkewCoordinates{:, 1}]));

int2Alt = 1/2*subs(int(diff(subs(wSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), eta) + subs(J22*thetaXSubs+J12*thetaYSubs, [xiBar; etaBar], [xi; eta] + (J0 \ c) * H1), eta), x1, xNodes(1));
deltaGamma22Alt = subs(int2Alt, [xi, eta], [1, -1]) - subs(int2Alt, [xi, eta], [-1, -1]);
deltaGamma23Alt = subs(int2Alt, [xi, eta], [1, 1]) - subs(int2Alt, [xi, eta], [-1, -1]);
deltaGamma24Alt = subs(int2Alt, [xi, eta], [-1, 1]) - subs(int2Alt, [xi, eta], [-1, -1]);
deltaGamma = M2*deltaGamma22+M3*deltaGamma23+M4*deltaGamma24;
GammaEtaBarDSG = diff(M2, etaBar)*deltaGamma22+diff(M3, etaBar)*deltaGamma23+diff(M4, etaBar)*deltaGamma24;
GammaEtaBarDSG = simplify(GammaEtaBarDSG);


GammaDSG = J0.'\[GammaXiBarDSG; GammaEtaBarDSG];

% Recreation:
GammaXiBarTest = subs(GammaXiBar, phi, 1);
GammaEtaBarTest = subs(GammaEtaBar, phi, 1);
% disp(GammaXiBarTest);
% disp(GammaEtaBarTest);

% There are multiple solutions to this equation!
% Therefore, plotting is necessary
% figure;
% fplot(subs(GammaXiBarTest, xiBar, 0));
% figure;
% fplot(subs(GammaEtaBarTest, etaBar, 0));

pointA = solve(subs(GammaXiBarTest, etaBar, 1), xiBar);
pointB = solve(subs(GammaEtaBarTest, xiBar, -1), etaBar);
pointC = solve(subs(GammaXiBarTest, etaBar, -1), xiBar);
pointD = solve(subs(GammaEtaBarTest, xiBar, 1), etaBar);

collocationPointsInSkewCoordinates = sym(zeros(2, 4));
collocationPointsInSkewCoordinates(:,1) = [simplify(pointA), 1];
collocationPointsInSkewCoordinates(:,2) = [-1, simplify(pointB)];
collocationPointsInSkewCoordinates(:,3) = [simplify(pointC), -1];
collocationPointsInSkewCoordinates(:,4) = [1, simplify(pointD)];
% Conclusion:
% The elimination of shear locking is not possible using the standard ANS
% interpolation. More precisely, to eliminate transverse shear locking in
% GammXiBar it is required to evaluate it at (xiBar, etaBar) = (0,
% arbitrary) and in GammaEtaBar at (xiBar, etaBar) = (arbitrary, -5).
% The next step is to check, whether these sampling points are mesh
% dependent (Assumption: yes).

% Results of further investigations:
% The elimination of transverse shear locking should be possible using a
% ANS-like interpolation. There is always a multitude of points, where the
% Kirchhoff condition is fulfilled. However, they are usually not lying on
% a line parallel to the xiBar/etaBar-Axis. Furthermore, these points are
% highly mesh dependent. Therefore, a quite complicated interpolation would
% be necessary.
% Question: Is this the same in the reference element?

% Further observation:
% The terms that cause transverse shear locking are just linear in xiBar
% and etaBar. This motivates an EAS approach.
