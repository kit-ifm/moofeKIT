% Script to investigate transverse shear locking for Bubnov-Galerkin FE

% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms w1 w2 w3 w4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi
syms xiBar etaBar

% plateLength
plateLength = 50;

% Nodal coordinates
xNodes = [0, 3/5*plateLength, 2/5*plateLength, 0];
yNodes = [0, 0, 2/3*plateLength, plateLength];

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

detJ = det(J);
detJ0 = det(J0);

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

wSubs = simplify(subs(w, [w1 w2 w3 w4], [0 (xNodes(2)^2/plateLength-xNodes(2))*phi (xNodes(3)^2/plateLength-xNodes(3))*phi 0]));

% Rotations
thetaX = M1*thetaX1 + M2*thetaX2 + M3*thetaX3 + M4*thetaX4;
thetaY = M1*thetaY1 + M2*thetaY2 + M3*thetaY3 + M4*thetaY4;

thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [phi (-2*xNodes(2)/plateLength+1)*phi (-2*xNodes(3)/plateLength+1)*phi phi]));


% Transverse shear strain
GammaXiBar = diff(wSubs, xiBar) + J011*thetaYSubs+J021*thetaXSubs;
GammaEtaBar = diff(wSubs, etaBar) + J022*thetaXSubs+J012*thetaYSubs;

% Backtransformation
Gamma = J0.'\[GammaXiBar; GammaEtaBar];


% Recreation:
GammaXiBarTest = subs(GammaXiBar, phi, 1);
GammaEtaBarTest = subs(GammaEtaBar, phi, 1);
disp(GammaXiBarTest);
disp(GammaEtaBarTest);

% There are multiple solutions to this equation!
% Therefore, plotting is necessary
figure;
fplot(subs(GammaXiBarTest, xiBar, 0));
figure;
fplot(subs(GammaEtaBarTest, etaBar, 0));

sol = solve(subs(GammaXiBarTest, xiBar, 0), etaBar);
disp(sol);

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
