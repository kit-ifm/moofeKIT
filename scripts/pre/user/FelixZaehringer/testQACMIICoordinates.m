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
% pointsInSkewCoordinates = [xi; eta] + [C1; C2] * H1;
nodalPointsInSkewCoordinates = sym(zeros(2, 4));
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates(1, ii) = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates(2, ii) = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

pointsInQACMIICoordinates = computeQACMIICoordinates([xi, eta], [xNodes; yNodes]);
nodalPointsInQACMIICoordinates = computeQACMIICoordinates(nodesXi, [xNodes; yNodes]);

% comparison
nodalPointsInQACMIICoordinatesEvalTest = eval(subs(nodalPointsInQACMIICoordinates, [xNodes; yNodes], nodesXi));
nodalPointsInQACMIICoordinatesEval = eval(subs(nodalPointsInQACMIICoordinates, [xNodes; yNodes], [0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));
nodalPointsInSkewCoordinatesEval = eval(subs(nodalPointsInSkewCoordinates, [xNodes; yNodes], [0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));

pointsInQACMIICoordinatesEval = eval(subs(pointsInQACMIICoordinates, [xNodes; yNodes], [0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));
pointsInSkewCoordinatesEval = eval(subs(pointsInSkewCoordinates, [xNodes; yNodes], [0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));

pointsInQACMIICoordinatesEvalTranslation = eval(subs(pointsInQACMIICoordinates, [xNodes; yNodes], [0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]+ones(2, 4)*23));
alpha=3.14;
pointsInQACMIICoordinatesEvalRotation = eval(subs(pointsInQACMIICoordinates, [xNodes; yNodes], [cosd(alpha), -sind(alpha); sind(alpha), cosd(alpha)]*[0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));
eval(simplify(pointsInQACMIICoordinatesEval - pointsInQACMIICoordinatesEvalRotation))

pointsInSkewCoordinatesEvalRotation = eval(subs(pointsInSkewCoordinates, [xNodes; yNodes], [cosd(alpha), -sind(alpha); sind(alpha), cosd(alpha)]*[0, 3, 2.7, 0.1; 0, 0.4, 3.4, 3.7]));
eval(simplify(pointsInSkewCoordinatesEval - pointsInSkewCoordinatesEvalRotation))

eval(simplify(pointsInQACMIICoordinatesEval - pointsInSkewCoordinatesEval))