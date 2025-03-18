syms xi eta;
syms xiBar etaBar;
syms C1 C2;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12;

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
nodalPointsInSkewCoordinates = cell(2, 4);
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates{1, ii} = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates{2, ii} = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end
pointsInSkewCoordinates = [xi; eta] + [C1; C2] * H1;
% pointsInSkewCoordinates = [xiBar; etaBar];


tau = [1, 0, pointsInSkewCoordinates(2), 0;...
       0, 1, 0, pointsInSkewCoordinates(1)];
% tau = [1, 0, eta, 0;...
%        0, 1, 0, xi];

% tau = adjoint(J)*J0*tau;

s = sym(zeros(size(tau)));
s(:, 1) = tau(:, 1);
s(:, 2) = tau(:, 2) - int(int(tau(:, 2).'*s(:, 1), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, 1).'*s(:, 1), xi, [-1, 1]), eta, [-1, 1]))*s(:, 1);

temp3 = sym(zeros(2, 1));
for k=1:2
temp3 = temp3 + int(int(tau(:, 3).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 3) = tau(:, 3) - temp3;

temp4 = sym(zeros(2, 1));
for k=1:3
temp4 = temp4 + int(int(tau(:, 4).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 4) = tau(:, 4) - temp4;

% check s
for i=1:4
    for j=1:4
        if j > i
            if i ~= j
                assert(int(int(subs(s(:, i).'*s(:, j), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);
            end
        end
    end
end

% compute gammaInit
w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;
GammaXiO = subs(diff(w, xi)+J011*thetaX+J021*thetaY, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1]);
GammaEtaO = subs(diff(w, eta)+J012*thetaX+J022*thetaY, [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1]);
[BsGammaXiO, BsGammaXiOT] = coeffs(GammaXiO, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);
[BsGammaEtaO, BsGammaEtaOT] = coeffs(GammaEtaO, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);

GammaXiANS = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
GammaEtaANS = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);
[BsGammaXiANS, BsGammaXiANST] = coeffs(GammaXiANS, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);
[BsGammaEtaANS, BsGammaEtaANST] = coeffs(GammaEtaANS, [w1, thetaX1, thetaY1, w2, thetaX2, thetaY2, w3, thetaX3, thetaY3, w4, thetaX4, thetaY4]);

GammaXiDiff = simplify(BsGammaXiANS-BsGammaXiO);
GammaEtaDiff = simplify(BsGammaEtaANS-BsGammaEtaO);

gammaInit = [GammaXiDiff; GammaEtaDiff];

% construct gamma
gammaInit = [xi, 0, xi*eta-1/3*C2*eta, 0; 0, eta, 0, xi*eta-1/3*C1*xi];
% gammaInit = J0.'*adjoint(J).'*gammaInit;

gamma = sym(zeros(size(gammaInit)));

for i=1:size(gammaInit, 2)
    tempGamma = sym(zeros(2, 1));
    for k=1:4
        tempGamma = tempGamma + int(int(gammaInit(:, i).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
    end
    gamma(:, i) = gammaInit(:, i) - tempGamma;
end

a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
h = 1/4*[1, -1, 1, -1].';
C1val = 1/detJ0*(a2.'*yNodes.'*h.'*xNodes.'-a2.'*xNodes.'*h.'*yNodes.');
C2val = 1/detJ0*(-a1.'*yNodes.'*h.'*xNodes.'+a1.'*xNodes.'*h.'*yNodes.');

% check gamma
for i=1:4
    for j=1:4
        assert(int(int(subs(subs(gammaInit(:, i).'*s(:, j), [C1, C2], [C1val, C2val]), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);  
    end
end

% Q Matrix
xVal = xNodes.';
yVal = yNodes.';
Q = sym(zeros(4, 12));
Q(1, [2, 11]) = a1.'*xVal;
Q(1, [3, 12]) = a1.'*yVal;
Q(1, [5, 8]) = -a1.'*xVal;
Q(1, [6, 9]) = -a1.'*yVal;
Q(2, [2, 5]) = a2.'*xVal;
Q(2, [3, 6]) = a2.'*yVal;
Q(2, [8, 11]) = -a2.'*xVal;
Q(2, [9, 12]) = -a2.'*yVal;
Q(3, [5, 11]) = a1.'*xVal;
Q(3, [6, 12]) = a1.'*yVal;
Q(3, [2, 8]) = -a1.'*xVal;
Q(3, [3, 9]) = -a1.'*yVal;
Q(4, [5, 11]) = a2.'*xVal;
Q(4, [6, 12]) = a2.'*yVal;
Q(4, [2, 8]) = -a2.'*xVal;
Q(4, [3, 9]) = -a2.'*yVal;
Q = 1/4*Q;

% special entries
M12 = -J0(2,1)*J0(2,2)-J0(1,1)*J0(1,2);
gammaInit = [xi-3*xi*eta^2, 0, xi*eta-1/3*C2*eta-1/3*C1*xi, 0; 0, eta-3*xi^2*eta, 0, xi*eta-1/3*C1*xi-1/3*C2*eta];

% check gamma transformed
GammaTransformed = detJ0/detJ*inv(J0.')*gammaInit*Q;
sTransformed = inv(J0.')*s;

for i=1:4
    for j=1:4
        assert(int(int(subs(subs(gammaInit(:, i).'*adjoint(J0)*sTransformed(:, j), [C1, C2], [C1val, C2val]), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);  
%         assert(int(int(subs(subs(GammaTransformed(:, i).'*sTransformed(:, j)*detJ, [C1, C2], [C1val, C2val]), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.2, 0.8, 1, -1, -1, -1, 1, 1]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);
%         assert(int(int(subs(subs(GammaTransformed(:, i).'*sTransformed(:, j)*detJ, [C1, C2], [C1val, C2val]), [x1, x2, x3, x4, y1, y2, y3, y4], [-1, 1, 1, -1, -1, -1, 1, 1]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);
%         assert(int(int(subs(subs(GammaTransformed(:, i).'*sTransformed(:, j)*detJ, [C1, C2], [C1val, C2val]), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);  
        disp(['ok: ', num2str(i), ', ', num2str(j)]);
    end
end


int(int(subs(gamma.'*tau, [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1]);


% compute Difference to ANS Interpolation
BsANSNew = [BsGammaXiO;BsGammaEtaO]+gamma;
BsANSDiff = simplify(BsANSNew - [BsGammaXiANS;BsGammaEtaANS])

% test
BsTest = BsANSNew - [BsGammaXiO;BsGammaEtaO];
int(int(subs(BsTest.'*tau, [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])
