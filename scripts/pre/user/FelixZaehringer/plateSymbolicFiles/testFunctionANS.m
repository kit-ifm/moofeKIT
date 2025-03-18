%% Quadrilateral
syms xi eta;
syms xiBar etaBar;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12;
A = sym('a', [8*6, 1]);
syms C1 C2;
syms a0x a1x a2x htx a0y a1y a2y hty;

a0 = 1/4*[1, 1, 1, 1];
a1 = 1/4*[-1, 1, 1, -1];
a2 = 1/4*[-1, -1, 1, 1];
h = 1/4*[1, -1, 1, -1];
xNodes = [x1, x2, x3, x4];
yNodes = [y1, y2, y3, y4];

N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);

x = a0x + a1x*xi + a2x*eta + htx*xi*eta;
y = a0y + a1y*xi + a2y*eta + hty*xi*eta;
% x = N1*x1+N2*x2+N3*x3+N4*x4;
% y = N1*y1+N2*y2+N3*y3+N4*y4;

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
pointsInSkewCoordinates = [xi; eta] + [C1; C2] * H1;


nodalPointsInSkewCoordinates = sym(zeros(2, 4));
for ii = 1:numberOfNodes
    nodalPointsInSkewCoordinates(1, ii) = simplify(subs(pointsInSkewCoordinates(1), {xi; eta}, nodesXi(:, ii)));
    nodalPointsInSkewCoordinates(2, ii) = simplify(subs(pointsInSkewCoordinates(2), {xi; eta}, nodesXi(:, ii)));
end

% pointsInSkewCoordinates = [xiBar; etaBar];


% Metric shape functions
continuumObject.elementGeometryType = 'quadrilateral';
[M_k_I, ~] = computeMetricShapeFunctions(continuumObject, 2, nodalPointsInSkewCoordinates, pointsInSkewCoordinates);
M1 = M_k_I(1);
M2 = M_k_I(2);
M3 = M_k_I(3);
M4 = M_k_I(4);

w = M1*w1+M2*w2+M3*w3+M4*w4;
thetaX = M1*thetaX1+M2*thetaX2+M3*thetaX3+M4*thetaX4;
thetaY = M1*thetaY1+M2*thetaY2+M3*thetaY3+M4*thetaY4;

GammaXiBarO = diff(w, xiBar)+J011*thetaX+J021*thetaY;
GammaEtaBarO = diff(w, etaBar)+J012*thetaX+J022*thetaY;

wIso = N1*w1+N2*w2+N3*w3+N4*w4;
thetaXIso = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaYIso = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(wIso, xi)+J11*thetaXIso+J21*thetaYIso;
GammaEtaO = diff(wIso, eta)+J12*thetaXIso+J22*thetaYIso;

[GammaXiOC, GammaXiOCT] = coeffs(GammaXiO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);
[GammaEtaOC, GammaEtaOCT] = coeffs(GammaEtaO, [w1 thetaX1 thetaY1 w2 thetaX2 thetaY2 w3 thetaX3 thetaY3 w4 thetaX4 thetaY4]);

BIso = [GammaXiOC; GammaEtaOC];

GammaXi_A = subs(GammaXiOC, [xi, eta], [0, 1]);
GammaXi_C = subs(GammaXiOC, [xi, eta], [0, -1]);
GammaEta_B = subs(GammaEtaOC, [xi, eta], [-1, 0]);
GammaEta_D = subs(GammaEtaOC, [xi, eta], [1, 0]);

BsANS = [1/2*(1+eta)*GammaXi_A+ 1/2*(1-eta)*GammaXi_C; 1/2*(1-xi)*GammaEta_B+ 1/2*(1+xi)*GammaEta_D];

S = [1, 0, pointsInSkewCoordinates(2), 0; 0, 1, 0, pointsInSkewCoordinates(1)];

% E = [1, 0, eta, 0; 0, 1, 0, xi];
% ENew = [xi^3, 0, xi, 0, eta^3, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2, 0; 0, xi^3, 0, eta^3, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2];
% 
% C = sym('c', [4, 4]);
% 
% equationOrthogonalityLeft = C.'*int(int(E.'*adjoint(J)*J0*S, xi, [-1, 1]), eta, [-1, 1]);
% equationOrthogonalityRightThetaX = int(int(BIso(:, 2:3:end).'*adjoint(J)*J0*S, xi, [-1, 1]), eta, [-1, 1]) - int(int(BsANS(:, 2:3:end).'*adjoint(J)*J0*S, xi, [-1, 1]), eta, [-1, 1]);
% equationOrthogonalityRightW= int(int(BIso(:, 1:3:end).'*adjoint(J)*J0*S, xi, [-1, 1]), eta, [-1, 1]) - int(int(BsANS(:, 1:3:end).'*adjoint(J)*J0*S, xi, [-1, 1]), eta, [-1, 1]);
% 
% 
% preFactorLocking = [1, -1, -1; 1, 1, -1; 1, 1, 1; 1, -1, 1];
% equationLocking = C*preFactorLocking;
% 
% equations = [equationOrthogonalityLeft(:) - equationOrthogonalityRight(:); equationLocking(:)];
% 
% sol = solve(equationOrthogonalityLeft(:) - equationOrthogonalityRightW(:), C);
% 
% % construct Ansatz Matrix
% F = sym(zeros(2, 4));
% numberFields = 6;
% for i=1:4
%     index1 = 1+(2*numberFields)*(i-1):numberFields+(2*numberFields)*(i-1);
%     index2 = numberFields+1+(2*numberFields)*(i-1):(2*numberFields)+(2*numberFields)*(i-1);
%     F(1, i) = [1, xi, eta, (1-eta^2), (1-xi^2), 1/2*(1-eta^2)*(1-xi^2)]*A(index1);
%     F(2, i) = [1, xi, eta, (1-eta^2), (1-xi^2), 1/2*(1-eta^2)*(1-xi^2)]*A(index2);
% %     F(1, i) = [1, xi, eta, xi*eta, xi^2, eta^2]*A(index1);
% %     F(2, i) = [1, xi, eta, xi*eta, xi^2, eta^2]*A(index2);
% %     F(1, i) = [1, eta, xi^2*eta^2]*A(index1);
% %     F(2, i) = [1, xi, xi^2*eta^2]*A(index2);
% end
% 
% preFactorMatrix = [GammaXi_A; GammaEta_B; GammaXi_C; GammaEta_D];
% % BsANS = F*preFactorMatrix;
% BsANS = [1/2*(1+eta)*GammaXi_A+ 1/2*(1-eta)*GammaXi_C; 1/2*(1-xi)*GammaEta_B+ 1/2*(1+xi)*GammaEta_D];
% 
% 
BsDiff = BIso - BsANS;
collect(BsDiff, [xi, eta])
% 
term = (BsDiff.'*adjoint(J)*J0).';
[termC, termT] = coeffs(term(1,12), [xi, eta]);

% BsDiffTOld = [1-eta^2, 0, xi, 0, xi*eta, 0, xi*eta^2, 0; 0, 1-xi^2, 0, eta, 0, xi*eta, 0, xi^2*eta];
% term = (BsDiffTOld.'*adjoint(J)*J0).';
% for i=1:8
%     [termC, termT] = coeffs(term(1,i), [xi, eta]);
%     disp(termT)
% end

% termT = [(eta^2-1), (xi^2-1), xi^2, xi, eta^2, eta, xi*eta, xi^2*eta, xi*eta^2, xi^2*eta^2];% evaluated from BsDiffOld
BsDiffT = sym(zeros(2, 16));
BsDiffT(1, 1:2:end) = flip(termT);
BsDiffT(2, 2:2:end) = flip(termT);

s = sym(zeros(size(BsDiffT)));
s(:, 1) = BsDiffT(:, 1);
s(:, 2) = BsDiffT(:, 2) - int(int(BsDiffT(:, 2).'*s(:, 1), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, 1).'*s(:, 1), xi, [-1, 1]), eta, [-1, 1]))*s(:, 1);

temp3 = sym(zeros(2, 1));
for k=1:2
temp3 = temp3 + int(int(BsDiffT(:, 3).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 3) = BsDiffT(:, 3) - temp3;


temp4 = sym(zeros(2, 1));
for k=1:3
temp4 = temp4 + int(int(BsDiffT(:, 4).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 4) = BsDiffT(:, 4) - temp4;

temp5 = sym(zeros(2, 1));
for k=1:4
temp5 = temp5 + int(int(BsDiffT(:, 5).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 5) = BsDiffT(:, 5) - temp5;

temp6 = sym(zeros(2, 1));
for k=1:5
temp6 = temp6 + int(int(BsDiffT(:, 6).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 6) = BsDiffT(:, 6) - temp6;

temp7 = sym(zeros(2, 1));
for k=1:6
temp7 = temp7 + int(int(BsDiffT(:, 7).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 7) = BsDiffT(:, 7) - temp7;

temp8 = sym(zeros(2, 1));
for k=1:7
temp8 = temp8 + int(int(BsDiffT(:, 8).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 8) = BsDiffT(:, 8) - temp8;

temp9 = sym(zeros(2, 1));
for k=1:8
temp9= temp9 + int(int(BsDiffT(:, 9).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 9) = BsDiffT(:, 9) - temp9;

temp10 = sym(zeros(2, 1));
for k=1:9
temp10= temp10 + int(int(BsDiffT(:, 10).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 10) = BsDiffT(:, 10) - temp10;

temp11 = sym(zeros(2, 1));
for k=1:10
temp11= temp11 + int(int(BsDiffT(:, 11).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 11) = BsDiffT(:, 11) - temp11;

temp12 = sym(zeros(2, 1));
for k=1:11
temp12= temp12 + int(int(BsDiffT(:, 12).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 12) = BsDiffT(:, 12) - temp12;

temp13 = sym(zeros(2, 1));
for k=1:12
temp13= temp13 + int(int(BsDiffT(:, 13).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 13) = BsDiffT(:, 13) - temp13;

temp14 = sym(zeros(2, 1));
for k=1:13
temp14= temp14 + int(int(BsDiffT(:, 14).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 14) = BsDiffT(:, 14) - temp14;

temp15 = sym(zeros(2, 1));
for k=1:14
temp15= temp15 + int(int(BsDiffT(:, 15).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 15) = BsDiffT(:, 15) - temp15;

temp16 = sym(zeros(2, 1));
for k=1:15
temp16= temp16 + int(int(BsDiffT(:, 16).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
s(:, 16) = BsDiffT(:, 16) - temp16;


% check s
for i=1:16
    for j=1:16
        if j > i
            if i ~= j
                assert(int(int(subs(s(:, i).'*s(:, j), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);
            end
        end
    end
end


% construct tau
tauInit = [1-3*xi^2-3*eta^2+9*xi^2*eta^2, 0, eta-5/3*eta^3, 0; 0, 1-3*xi^2-3*eta^2+9*xi^2*eta^2, 0, xi-5/3*xi^3];

tau = sym(zeros(size(tauInit)));

tempGamma1 = sym(zeros(2, 1));
for k=1:16
    tempGamma1 = tempGamma1 + int(int(tauInit(:, 1).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
tau(:, 1) = tauInit(:, 1) - tempGamma1;
disp(tau);

tempGamma2 = sym(zeros(2, 1));
for k=1:16
    tempGamma2 = tempGamma2 + int(int(tauInit(:, 2).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
tau(:, 2) = tauInit(:, 2) - tempGamma2;

tempGamma3 = sym(zeros(2, 1));
for k=1:16
    tempGamma3 = tempGamma3 + int(int(tauInit(:, 3).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
tau(:, 3) = tauInit(:, 3) - tempGamma3;

tempGamma4 = sym(zeros(2, 1));
for k=1:16
    tempGamma4 = tempGamma4 + int(int(tauInit(:, 4).'*s(:, k), xi, [-1, 1]), eta, [-1, 1])/(int(int(s(:, k).'*s(:, k), xi, [-1, 1]), eta, [-1, 1]))*s(:, k);
end
tau(:, 4) = tauInit(:, 4) - tempGamma4;


% check gamma
for i=1:4
    for j=1:16
        assert(int(int(subs(tau(:, i).'*s(:, j), [x1, x2, x3, x4, y1, y2, y3, y4], [-1.5, 1, 1.3, -1.1, -1.5, -0.9, 1.1, 1.2]), xi, [-1, 1]), eta, [-1, 1])==0, ['error for i=', num2str(i), ', j=', num2str(j), '.']);  
    end
end






