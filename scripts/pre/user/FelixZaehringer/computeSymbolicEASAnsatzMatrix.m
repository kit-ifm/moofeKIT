% Computes analytically the EAS Ansatz Matrix for plate elements

clear;
clc;

syms xi  eta;
syms x1 x2 x3 x4 y1 y2 y3 y4;
syms w1 w2 w3 w4 bx1 bx2 bx3 bx4 by1 by2 by3 by4;

nodalPoints = [x1, x2, x3, x4; y1, y2, y3, y4];
isoparametricNodes = [-1, 1, 1, -1; -1, -1, 1, 1];
arbitraryNodes = [0, 3, 4, 2; -1, -1, 3, 1];
u = [w1, bx1, by1, w2, bx2, by2, w3, bx3, by3, w4, bx4, by4].';
w = [w1, w2, w3, w4].';
thetaX = [bx1, bx2, bx3, bx4].';
thetaY = [by1, by2, by3, by4].';

% Lagrange shape functions
numberOfNodes = 4;
NSymbolic(1) = (1 - xi) .* (1 - eta) / 4;
NSymbolic(2) = (1 + xi) .* (1 - eta) / 4;
NSymbolic(3) = (1 + xi) .* (1 + eta) / 4;
NSymbolic(4) = (1 - xi) .* (1 + eta) / 4;

dNSymbolic = cell(2, 4);
for ii=1:numberOfNodes
    dNSymbolic{1, ii} = diff(NSymbolic(ii), xi);
    dNSymbolic{2, ii} = diff(NSymbolic(ii), eta);
end

% geometry coordinates
x=0;
y=0;
for ii=1:numberOfNodes
    x = x + NSymbolic(ii)*nodalPoints(1, ii);
    y = y + NSymbolic(ii)*nodalPoints(2, ii);
end

% geometry derivatives
x_xi = diff(x, xi);
x_eta = diff(x, eta);
y_xi = diff(y, xi);
y_eta = diff(y, eta);

% Jacobian
J = [x_xi, x_eta; y_xi, y_eta];
J0 = subs(J, {xi, eta}, {0, 0});

% simplification of Jacobian
syms J0_11 J0_12 J0_21 J0_22 hTx hTy;
J = [J0_11+eta*hTx, J0_12+xi*hTx; J0_21+eta*hTy, J0_22+xi*hTy];
J0 = subs(J, {xi, eta}, {0, 0});

detJ = det(J);
detJ0 = det(J0);

% compute EAS Ansatz
% syms ErU10 ErU11 ErU12 ErU13 ErU30 ErU31 ErU32 ErU33 ErU60 ErU61 ErU62 ErU63 ErU80 ErU81 ErU82 ErU83;
% syms ErU10 ErU11 ErU12 ErU13 ErU33 ErU62 ErU83;
% unknownVars = [ErU10 ErU11 ErU12 ErU13 ErU33 ErU62 ErU83];
% syms ErU10 ErU11 ErU12 ErU13 ErU30 ErU31 ErU32 ErU33 ErU60 ErU61 ErU62 ErU63 ErU80 ErU81 ErU82 ErU83;
% unknownVars = [ErU10 ErU11 ErU12 ErU13 ErU30 ErU31 ErU32 ErU33 ErU60 ErU61 ErU62 ErU63 ErU80 ErU81 ErU82 ErU83];
% Er = [ErU10 + ErU11*xi + ErU12*eta + ErU13*xi*eta, 0, ErU30 + ErU31*xi + ErU32*eta + ErU33*xi*eta, 0; 0, ErU60 + ErU61*xi + ErU62*eta + ErU63*xi*eta, 0, ErU80 + ErU81*xi + ErU82*eta + ErU83*xi*eta];
% Er = [ErU10 + ErU11*xi + ErU12*eta + ErU13*xi*eta, 0, ErU33*xi*eta, 0; 0, ErU62*eta, 0, ErU83*xi*eta];
Er = [xi, 0, xi*eta, 0, xi*eta^2, 0; 0, eta, 0, xi*eta, 0, xi^2*eta];
% syms ErU1 ErU2 ErU3 ErU4;
% Er = [ErU1, 0, ErU2, 0; 0, ErU3, 0, ErU4];
G = detJ0/detJ*J0.'\Er;

% syms ErU30 ErU31 ErU32 ErU33 ErU80 ErU81 ErU82 ErU83;
% unknownVars = [ErU30 ErU31 ErU32 ErU33 ErU80 ErU81 ErU82 ErU83];
% ErPetrovGalerkin = [ErU10 + ErU11*xi + ErU12*eta + ErU13*xi*eta, 0, ErU30 + ErU31*xi + ErU32*eta + ErU33*xi*eta, 0; 0, ErU60 + ErU61*xi + ErU62*eta + ErU63*xi*eta, 0, ErU80 + ErU81*xi + ErU82*eta + ErU83*xi*eta];
% ErPetrovGalerkin = [xi, 0, ErU30 + ErU31*xi + ErU32*eta + ErU33*xi*eta, 0, xi*eta^2, 0; 0, eta, 0, ErU80 + ErU81*xi + ErU82*eta + ErU83*xi*eta, 0, xi^2*eta];
% ErPetrovGalerkin = [pointsInSkewCoordinates(1), 0, ErU1, 0; 0, pointsInSkewCoordinates(2), 0, ErU2];
% H = detJ0/detJ*J0.'\ErPetrovGalerkin; % ist diese Transformation so erlaubt?!

% compute GammaBar
a0 = 1/4*[1, 1, 1, 1].';
a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
a3 = 1/4*[1, -1, 1, -1].';
% omega1 = a1.'*w;
% omega2 = a2.'*w;
% omega3 = a3.'*w;
% alpha0 = a0.'*thetaX;
% alpha1 = a1.'*thetaX;
% alpha2 = a2.'*thetaX;
% alpha3 = a3.'*thetaX;
% beta0 = a0.'*thetaY;
% beta1 = a1.'*thetaY;
% beta2 = a2.'*thetaY;
% beta3 = a3.'*thetaY;

syms omega1 omega2 omega3 alpha0 alpha1 alpha2 alpha3 beta0 beta1 beta2 beta3;

xNodes = nodalPoints(1, :).';
yNodes = nodalPoints(2, :).';

% GammaBar(1) = (J0(1,1)+a3.'*nodalPoints(1, :).'*eta)*(alpha0+alpha1*xi+alpha2*eta+alpha3*xi*eta)+(J0(2,1)+a3.'*nodalPoints(2, :).'*eta)*(beta0+beta1*xi+beta2*eta+beta3*xi*eta)+omega1+omega3*eta;
% GammaBar(2) = (J0(1,2)+a3.'*nodalPoints(1, :).'*xi)*(alpha0+alpha1*xi+alpha2*eta+alpha3*xi*eta)+(J0(2,2)+a3.'*nodalPoints(2, :).'*xi)*(beta0+beta1*xi+beta2*eta+beta3*xi*eta)+omega2+omega3*xi;
% hTx = a3.'*xNodes;
% hTy = a3.'*yNodes;

GammaBar(1) = (J0(1,1)+hTx*eta)*(alpha0+alpha1*xi+alpha2*eta+alpha3*xi*eta)+(J0(2,1)+hTy*eta)*(beta0+beta1*xi+beta2*eta+beta3*xi*eta)+omega1+omega3*eta;
GammaBar(2) = (J0(1,2)+hTx*xi)*(alpha0+alpha1*xi+alpha2*eta+alpha3*xi*eta)+(J0(2,2)+hTy*xi)*(beta0+beta1*xi+beta2*eta+beta3*xi*eta)+omega2+omega3*xi;

% ANS interpolation
GammaBarCollA = subs(GammaBar, {xi, eta}, {0, 1});
GammaBarCollB = subs(GammaBar, {xi, eta}, {-1, 0});
GammaBarCollC = subs(GammaBar, {xi, eta}, {0, -1});
GammaBarCollD = subs(GammaBar, {xi, eta}, {1, 0});

GammaBarANS(1) = collect(simplify(collect(1/2*((1+eta)*GammaBarCollA(1)+(1-eta)*GammaBarCollC(1)), [xi, eta])), [xi, eta]);
GammaBarANS(2) = collect(simplify(collect(1/2*((1+xi)*GammaBarCollD(2)+(1-xi)*GammaBarCollB(2)), [xi, eta])), [xi, eta]);

syms omegaMetric1 omegaMetric2 omegaMetric3 alphaMetric0 alphaMetric1 alphaMetric2 alphaMetric3 betaMetric0 betaMetric1 betaMetric2 betaMetric3;

metricCoord = [xi; eta] + inv(J0)*1/4*[hTx; hTy]*xi*eta;
xiBar = metricCoord(1);
etaBar = metricCoord(2);

GammaBarMetric(1) = J0(1,1)*(alphaMetric0+alphaMetric1*xiBar+alphaMetric2*etaBar+alphaMetric3*xiBar*etaBar)+J0(2,1)*(betaMetric0+betaMetric1*xiBar+betaMetric2*etaBar+betaMetric3*xiBar*etaBar)+omegaMetric1+omegaMetric3*etaBar;
GammaBarMetric(2) = J0(1,2)*(alphaMetric0+alphaMetric1*xiBar+alphaMetric2*etaBar+alphaMetric3*xiBar*etaBar)+J0(2,2)*(betaMetric0+betaMetric1*xiBar+betaMetric2*etaBar+betaMetric3*xiBar*etaBar)+omegaMetric2+omegaMetric3*xiBar;


% compute Bs
dN_X_I = J.' \ dNSymbolic;

Bs(1, 1:3:12) = dN_X_I(1, :);
Bs(2, 1:3:12) = dN_X_I(2, :);
Bs(1, 2:3:12) = NSymbolic;
Bs(2, 3:3:12) = NSymbolic;

test = simplify(Bs*u - J.'\GammaBar.', 'steps', 20); % passt

polynomials = collect(collect([J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBar.', [beta0, beta1, beta2, beta3]), [alpha0, alpha1, alpha2, alpha3]);

%solve for ANS modes
syms EA10 EA11 EA12 EA13 EA14 EA15 EA16 EA17 EA20 EA21 EA22 EA23 EA24 EA25 EA26 EA27;
unknowns = [EA10 EA20 EA11 EA21 EA12 EA22 EA13 EA23 EA14 EA24 EA15 EA25 EA16 EA26 EA17 EA27];
EA(1) = EA10+EA11*xi+EA12*eta+EA13*xi*eta+EA14*xi^2+EA15*eta^2+EA16*xi^2*eta+EA17*xi*eta^2;
EA(2) = EA20+EA21*xi+EA22*eta+EA23*xi*eta+EA24*xi^2+EA25*eta^2+EA26*xi^2*eta+EA27*xi*eta^2;
T1 = collect([J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBar.', [xi, eta]);
T2 = collect([J0(2,2), -J0(2,1); -J0(1,2), J0(1,1)]*EA.', [xi, eta]);
T3 = collect([J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBarANS.', [xi, eta]);
T = collect(T1+T2-T3, [xi, eta]);
[coeff1, terms1] = coeffs(T(1), [xi, eta]);
[coeff2, terms2] = coeffs(T(2), [xi, eta]);
sol = solve([coeff1; coeff2], unknowns);

% solve for ANS modes (metric shape functions)
syms EA18 EA28 EA19 EA29 EA110 EA210 EA111 EA211 EA112 EA212;
unknowns = [EA10 EA20 EA11 EA21 EA12 EA22 EA13 EA23 EA14 EA24 EA15 EA25 EA16 EA26 EA17 EA27 EA18 EA28 EA19 EA29 EA110 EA210 EA111 EA211 EA112 EA212];
EA(1) = EA10+EA11*xi+EA12*eta+EA13*xi*eta+EA14*xi^2+EA15*eta^2+EA16*xi^2*eta+EA17*xi*eta^2+EA18*xi^2*eta^2+EA19*xi^3*eta+EA110*xi*eta^3+EA111*xi^3*eta^2+EA112*xi^2*eta^3;
EA(2) = EA20+EA21*xi+EA22*eta+EA23*xi*eta+EA24*xi^2+EA25*eta^2+EA26*xi^2*eta+EA27*xi*eta^2+EA28*xi^2*eta^2+EA29*xi^3*eta+EA210*xi*eta^3+EA211*xi^3*eta^2+EA212*xi^2*eta^3;

% unknowns = [EA10 EA20 EA11 EA21 EA12 EA22 EA13 EA23 EA14 EA24 EA15 EA25 EA16 EA26 EA17 EA27];
% EA(1) = EA10+EA11*xiBar+EA12*etaBar+EA13*xiBar*etaBar+EA14*xiBar^2+EA15*etaBar^2+EA16*xiBar^2*etaBar+EA17*xiBar*etaBar^2;
% EA(2) = EA20+EA21*xiBar+EA22*etaBar+EA23*xiBar*etaBar+EA24*xiBar^2+EA25*etaBar^2+EA26*xiBar^2*etaBar+EA27*xiBar*etaBar^2;

T1 = collect(detJ*inv(J0.')*GammaBarMetric.', [xi, eta]);
T2 = collect([J0(2,2), -J0(2,1); -J0(1,2), J0(1,1)]*EA.', [xi, eta]);
T3 = collect([J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBarANS.', [xi, eta]);
T = collect(T1+T2-T3, [xi, eta]);
[coeff1, terms1] = coeffs(T(1), [xi, eta]);
[coeff2, terms2] = coeffs(T(2), [xi, eta]);
sol = solve([coeff1; coeff2], unknowns);

% compute E and alpha
E = [1, 0, xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2, 0, xi^2*eta^2, 0, xi^3*eta, 0, xi*eta^3, 0, xi^3*eta^2, 0, xi^2*eta^3, 0; 0, 1, 0, xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2, 0, xi^2*eta^2, 0, xi^3*eta, 0, xi*eta^3, 0, xi^3*eta^2, 0, xi^2*eta^3];



% other solutions (Lagrange shape functions)
prePolyn = collect(J0.'*polynomials, [xi, eta]);

% E = [xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2, 0; 0, xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2];
% 
% alpha = -[J0(1,1)*alpha1;
%     J0(1,2)*alpha1;
%     J0(2,1)*beta2;
%     J0(2,2)*beta2;
%     1/detJ0*(J0(1,1)*((J0(2,2)*hTx-J0(1,2)*hTy)*alpha1+alpha3*detJ0)+J0(2,1)*((J0(1,1)*hTy-J0(2,1)*hTx)*beta2+beta3*detJ0));
%     1/detJ0*(J0(1,2)*((J0(2,2)*hTx-J0(1,2)*hTy)*alpha1+alpha3*detJ0)+J0(2,2)*((J0(1,1)*hTy-J0(2,1)*hTx)*beta2+beta3*detJ0));
%     J0(1,1)/detJ0*(J0(1,1)*hTy-J0(2,1)*hTx)*alpha1;
%     J0(1,2)/detJ0*(J0(1,1)*hTy-J0(2,1)*hTx)*alpha1;
%     J0(2,1)/detJ0*(J0(2,2)*hTx-J0(1,2)*hTy)*beta2;
%     J0(2,2)/detJ0*(J0(2,2)*hTx-J0(1,2)*hTy)*beta2;
%     1/detJ0*(J0(1,1)*hTy-J0(2,1)*hTx)*(J0(1,1)*alpha3+J0(2,1)*beta3);
%     1/detJ0*(J0(1,1)*hTy-J0(2,1)*hTx)*(J0(1,2)*alpha3+J0(2,2)*beta3);
%     1/detJ0*(J0(2,2)*hTx-J0(1,2)*hTy)*(J0(1,1)*alpha3+J0(2,1)*beta3);
%     1/detJ0*(J0(2,2)*hTx-J0(1,2)*hTy)*(J0(1,2)*alpha3+J0(2,2)*beta3);
%     ];

E = [1, 0, xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2, 0; 0, 1, 0, xi, 0, eta, 0, xi*eta, 0, xi^2, 0, eta^2, 0, xi^2*eta, 0, xi*eta^2];

alpha = [sol.EA10;
    sol.EA20;
    sol.EA11;
    sol.EA21;
    sol.EA12;
    sol.EA22;
    sol.EA13;
    sol.EA23;
    sol.EA14;
    sol.EA24;
    sol.EA15;
    sol.EA25;
    sol.EA16;
    sol.EA26;
    sol.EA17;
    sol.EA27];

test = simplify([J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBar.'+[J0(2,2), -J0(2,1); -J0(1,2), J0(1,1)]*E*alpha-[J(2,2), -J(2,1); -J(1,2), J(1,1)]*GammaBarANS.', 'steps', 20);
test2 = simplify(J.'\GammaBar.' + detJ0/detJ*inv(J0.')*E*alpha-J.'\GammaBarANS.', 'steps', 20);

% epsLocal = collect(simplify(collect(GammaBar.' + E*alpha, [xi, eta])), [xi, eta]);
eps = collect(simplify(collect(collect(simplify(J.'\GammaBar.' + detJ0/detJ*inv(J0.')*E*alpha), [alpha0, alpha1, alpha2, alpha3]), [beta0, beta1, beta2, beta3])), [xi, eta]);

Eps = collect(J.'*eps, [xi, eta]);
test3 = simplify(Eps - GammaBarANS.');
% eps = collect(simplify(collect(simplify(eps), [xi, eta])), [xi, eta]);

epsTrueSol = [(J0(2,2)*J0(1,1)+J0(2,2)*hTx*eta+J0(1,1)*hTy*xi-J0(2,1)*J0(1,2)-J0(2,1)*hTx*xi-J0(1,2)+hTy*eta)*xi*alpha1+(J0(2,2)*J0(1,1)+J0(2,2)*hTx*eta+J0(1,1)*hTy*xi-J0(2,1)*J0(1,2)-J0(2,1)*hTx*xi-J0(1,2)+hTy*eta)*xi*eta*alpha3];
test1 = collect(eps(1)-epsTrueSol(1), [alpha0, alpha1, alpha2, alpha3]);
test2 = collect(eps(2)-epsTrueSol(2), [beta0, beta1, beta2, beta3]);
% compute right side
disp('Right Side...');
% detJ = subs(detJ, nodalPoints, isoparametricNodes);
% G = subs(G, nodalPoints, isoparametricNodes);
% H = subs(H, nodalPoints, isoparametricNodes);
% Bs = subs(Bs, nodalPoints, isoparametricNodes);
disp('- T1');
T1 = collect(simplify(G.'*Bs*u*detJ), [xi, eta]);
disp('- T2');
T2 = collect(simplify(G.'*H*detJ), [xi, eta]);
disp('- T1Int');
T1Int = int(int(T1, xi, -1, 1), eta, -1, 1);
disp('- T2Int');
T2Int = int(int(T2, xi, -1, 1), eta, -1, 1);
% T1Int = int(int(T1, xi), eta);
% T2Int = int(int(T2, xi), eta);
% disp('- T3');
% T3 = T2Int \ T1Int;
% rightSide = -T3;

% testRightSide = subs(rightSide, nodalPoints, isoparametricNodes);

save EASRightSide.mat rightSide;

% compute left side
disp('Left Side...');
a0 = 1/4*[1, 1, 1, 1].';
a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
a3 = 1/4*[1, -1, 1, -1].';
alpha1 = a1.'*thetaX;
alpha2 = a2.'*thetaX;
alpha3 = a3.'*thetaX;
beta1 = a1.'*thetaY;
beta2 = a2.'*thetaY;
beta3 = a3.'*thetaY;
leftSide(1) = J0_11*alpha1+J0_21*beta1;
leftSide(2) = J0_12*alpha2+J0_22*beta2;
leftSide(3) = J0_11*alpha3+J0_21*beta3+hTx*alpha1+hTy*beta1;
leftSide(4) = J0_12*alpha3+J0_22*beta3+hTx*alpha2+hTy*beta2;
leftSide(5) = hTx*alpha3+hTy*beta3;
leftSide(6) = hTx*alpha3+hTy*beta3;
leftSide = leftSide.';
% testLeftSide = subs(leftSide, nodalPoints, isoparametricNodes);
% 
disp('Equation...')
equation = collect(simplify(T2Int*leftSide-T1Int), unknownVars);

save EASEquation.mat T1Int T2Int leftSide equation;

disp('Solving...');
% % subs(equation, {ErU10 ErU11 ErU12 ErU13 ErU30 ErU31 ErU32 ErU33 ErU60 ErU61 ErU62 ErU63 ErU80 ErU81 ErU82 ErU83}, {0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 1});
% % subs(equation, unknownVars, {0 0 0 0 0 0 0});
% % S = vpasolve(equation, [ErU10 ErU11 ErU12 ErU13 ErU30 ErU31 ErU32 ErU33 ErU60 ErU61 ErU62 ErU63 ErU80 ErU81 ErU82 ErU83], [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 1; 0; 0; 0; 0; 1]);
S = solve(equation, unknownVars);
disp('Done.');
% % test for isoparametric element
% testG = subs(G, nodalPoints, isoparametricNodes);