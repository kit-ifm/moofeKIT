% Computes analytically the EAS Ansatz Matrix for plate elements

clear;
clc;

syms xi  eta;
syms x1 x2 x3 x4 y1 y2 y3 y4;

nodalPoints = {x1, x2, x3, x4; y1, y2, y3, y4};
isoparametricNodes = [-1, 1, 1, -1; -1, -1, 1, 1];

% Lagrange shape functions
numberOfNodes = 4;
NSymbolic{1} = (1 - xi) .* (1 - eta) / 4;
NSymbolic{2} = (1 + xi) .* (1 - eta) / 4;
NSymbolic{3} = (1 + xi) .* (1 + eta) / 4;
NSymbolic{4} = (1 - xi) .* (1 + eta) / 4;

dNSymbolic = cell(2, 4);
for ii=1:numberOfNodes
    dNSymbolic{1, ii} = diff(NSymbolic{ii}, xi);
    dNSymbolic{2, ii} = diff(NSymbolic{ii}, eta);
end

% geometry coordinates
x=0;
y=0;
for ii=1:numberOfNodes
    x = x + NSymbolic{ii}*nodalPoints{1, ii};
    y = y + NSymbolic{ii}*nodalPoints{2, ii};
end

% geometry derivatives
x_xi = diff(x, xi);
x_eta = diff(x, eta);
y_xi = diff(y, xi);
y_eta = diff(y, eta);

% Jacobian
J = [x_xi, x_eta; y_xi, y_eta];
J0 = subs(J, {xi, eta}, {0, 0});

detJ = det(J);
detJ0 = det(J0);

% compute skew coordinates
nodesXi = elementNodesInLocalCoordinates(2, 'quadrilateral', 4);
[hourglassVector, ~] = computeHourglassVectorAndFunction(2, [xi; eta]);
c = {0; 0};
for ii = 1:numberOfNodes
    c{1} = c{1} + 1 / numberOfNodes * nodalPoints{1, ii} * hourglassVector(ii);
    c{2} = c{2} + 1 / numberOfNodes * nodalPoints{2, ii} * hourglassVector(ii);
end
H1 = xi * eta;
pointsInSkewCoordinates = [xi; eta] + (J0 \ c) * H1;
nodalPointsInSkewCoordinates = cell(2, numberOfNodes);
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates{1, ii} = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates{2, ii} = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

% compute metric shape functions
xPoints = [nodalPointsInSkewCoordinates{1, :}].';
yPoints = [nodalPointsInSkewCoordinates{2, :}].';
P = [ones(4, 1), xPoints, yPoints, xPoints .* yPoints];
p = [1, pointsInSkewCoordinates(1), pointsInSkewCoordinates(2), pointsInSkewCoordinates(1) .* pointsInSkewCoordinates(2)];
M_k_I = p / P;

for ii=1:numberOfNodes
   dMr(1, ii) = diff(M_k_I(ii), xi);
   dMr(2, ii) = diff(M_k_I(ii), eta);
end

% compute EAS Ansatz
Er = [xi, 0, xi*eta, 0; 0, eta, 0, xi*eta];
G = detJ0/detJ*J.'\Er;

syms ErU1 ErU2;
ErPetrovGalerkin = [pointsInSkewCoordinates(1), 0, ErU1, 0; 0, pointsInSkewCoordinates(2), 0, ErU2];
H = detJ0/detJ*J.'\ErPetrovGalerkin; % ist diese Transformation so erlaubt?!


% compute Ls
dM_X_I = J0.' \ dMr;

Ls(1, 1:3:12) = dM_X_I(1, :);
Ls(2, 1:3:12) = dM_X_I(2, :);
Ls(1, 2:3:12) = M_k_I;
Ls(2, 3:3:12) = M_k_I;

% compute right side
rightSide = int(G.'*Ls*detJ, xi, [-1, 1]);

% test for isoparametric element
testG = subs(G, nodalPoints, isoparametricNodes);
testH = subs(H, nodalPoints, isoparametricNodes);