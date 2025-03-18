% Script to investigate transverse shear locking for Petrov-Galerkin FE

% Syms
syms xi eta
syms xiBar etaBar
syms x1 y1 x2 y2 x3 y3 x4 y4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi

% bendingType: bending around which axis
bendingType = 'y';

% elementForm
distortionType = 'none';

% Nodal coordinates
if strcmp(distortionType, 'none')
    xNodes = [-1, 1, 1, -1];
    yNodes = [-1, -1, 1, 1];
elseif strcmp(distortionType, 'parallel')
    if strcmp(bendingType, 'x')
        xNodes = [-1, 1, 2, 0];
        yNodes = [-1, -1, 1, 1];
    elseif strcmp(bendingType, 'y')
        xNodes = [-1, 1, 1, -1];
        yNodes = [-1, -2, 0, 1];
    end
elseif strcmp(distortionType, 'angular')
    if strcmp(bendingType, 'x')
        xNodes = [-2, 3, 1, -1];
        yNodes = [-1, -1, 1, 1];
    elseif strcmp(bendingType, 'y')
        xNodes = [-1, 1, 1, -1];
        yNodes = [-2, -1, 1, 3];
    end
end
numberOfNodes = 4;

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

% Jacobian
J = [J11, J12; J21, J22];
J0 = subs(J, {xi, eta}, {0, 0});

detJ = det(J);
detJ0 = det(J0);

% skew coordinates
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

% Rotations
thetaX = M1*thetaX1 + M2*thetaX2 + M3*thetaX3 + M4*thetaX4;
thetaY = M1*thetaY1 + M2*thetaY2 + M3*thetaY3 + M4*thetaY4;

if strcmp(bendingType, 'x')
    thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [phi phi -phi -phi]));
    thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [0 0 0 0]));
elseif strcmp(bendingType, 'y')
    thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
    thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [phi -phi -phi phi]));
end

% Transverse shear strain
GammaX = J0(1,1)*thetaYSubs+J0(2,1)*thetaXSubs;
GammaY = J0(2,2)*thetaXSubs+J0(1,2)*thetaYSubs;

% Find matching xi/eta for Kirchhoff hypothesis
disp('Start solving...');
sol = solve([GammaX, GammaY], [xiBar, etaBar]);
disp(GammaX);
disp(GammaY);
disp(sol);

% Conclusion:
% The Kirchhoff condition is satisfied for any element geometry, when
% interpolating GammaX at xi=0 (eta is arbitrary) and GammaY at eta=0 (xi
% is arbitrary).

% Correct solution
GammaX = subs(J0(1,1)*thetaY+J0(2,1)*thetaX, xiBar, 0);
GammaY = subs(J0(2,2)*thetaX+J0(1,2)*thetaY, etaBar, 0);
disp(GammaX);
disp(GammaY);

% ANS solution
GammaXANS = 1/2*(1+etaBar)*subs(J0(1,1)*thetaY+J0(2,1)*thetaX, [xiBar, etaBar], [0, 1])+1/2*(1-etaBar)*subs(J0(1,1)*thetaY+J0(2,1)*thetaX, [xiBar, etaBar], [0, -1]);

disp(GammaXANS);
disp(simplify(GammaX-GammaXANS))
