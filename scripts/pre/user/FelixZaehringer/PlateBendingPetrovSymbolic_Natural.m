% Script to investigate transverse shear locking for Petrov-Galerkin FE

% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms w1 w2 w3 w4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi

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

% verification point: Jacobian => ok
J_check = eval(subs(J, [xi, eta], [-0.5774, -0.5774]));

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

% verification point: dXiBar_xi => ok
dXiBar_xi = diff(pointsInSkewCoordinates, xi);
dXiBar_eta = diff(pointsInSkewCoordinates, eta);
dXiBar_Xi = eval(subs([dXiBar_xi, dXiBar_eta], [xi, eta], [-0.5774, -0.5774]));

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

% verification point: metric shape functions => ok
M_check_A = eval(subs([M1 M2 M3 M4], [xi, eta], [0, 1]));
M_check_B = eval(subs([M1 M2 M3 M4], [xi, eta], [-1, 0]));
M_check_C = eval(subs([M1 M2 M3 M4], [xi, eta], [0, -1]));
M_check_D = eval(subs([M1 M2 M3 M4], [xi, eta], [1, 0]));

% verification point: derivatives
dM_check_A = eval(subs(diff([M1 M2 M3 M4], xi), [xi, eta], [0, 1]));

% displacement
w = M1*w1 + M2*w2 + M3*w3 + M4*w4;

wSubs = simplify(subs(w, [w1 w2 w3 w4], [0 (xNodes(2)^2/plateLength-xNodes(2))*phi (xNodes(3)^2/plateLength-xNodes(3))*phi 0]));

% Rotations
thetaX = M1*thetaX1 + M2*thetaX2 + M3*thetaX3 + M4*thetaX4;
thetaY = M1*thetaY1 + M2*thetaY2 + M3*thetaY3 + M4*thetaY4;

thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [phi (-2*xNodes(2)/plateLength+1)*phi (-2*xNodes(3)/plateLength+1)*phi phi]));


% Transverse shear strain
GammaXi = diff(wSubs, xi) + J11*thetaYSubs+J21*thetaXSubs;
GammaEta = diff(wSubs, eta) + J22*thetaXSubs+J12*thetaYSubs;

% ANS solution
GammaXiANS = 1/2*(1+eta)*subs(GammaXi, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXi, [xi, eta], [0, -1]);
GammaEtaANS = 1/2*(1+xi)*subs(GammaEta, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEta, [xi, eta], [-1, 0]);

% Backtransformation
Gamma = J.'\[GammaXi; GammaEta];
GammaANS = J.'\[GammaXiANS; GammaEtaANS];

% Recreation:
GammaXiTest = subs(GammaXi, phi, 1);
GammaEtaTest = subs(GammaEta, phi, 1);
disp(GammaXiTest);
disp(GammaEtaTest);

% There are multiple solutions to this equation!
% Therefore, plotting is necessary
figure;
fplot(subs(GammaXiTest, xi, 0));
figure;
fsurf(GammaXiTest);
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])
figure;
fplot(subs(GammaEtaTest, eta, 0));
% 
sol = solve(subs(GammaXiTest, xi, 0), eta);
disp(sol);

% Conclusion:
% The elimination of shear locking seems to be possible using the standard ANS
% interpolation. More precisely, to eliminate transverse shear locking in
% GammXi it is required to evaluate it at (xi, eta) = (0, 1/-1) and in
% GammaEta at (xi, eta) = (1/-1, 0).
% This needs to be investigated more.

pointsInSkewCoordinatesEval1 = subs(pointsInSkewCoordinates, [xi, eta], [0, 1]);
pointsInSkewCoordinatesEval2 = subs(pointsInSkewCoordinates, [xi, eta], [0, -1]);

GammaXiBarTest = J0.'*inv(J.')*[subs(GammaXiTest, [xi, eta], [0, 1]); subs(GammaEtaTest, [xi, eta], [0, 1])];

figure;
for xiVar=-1:0.1:1
    for etaVar = -1:0.1:1
        if subs(subs(GammaXi, phi, 1), [xi, eta], [xiVar, etaVar]) == 0
            disp(xiVar);
            disp(etaVar);
        end
    end
end