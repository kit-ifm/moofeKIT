% Script to investigate transverse shear locking for Bubnov-Galerkin FE

% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms w1 w2 w3 w4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi omega MEI
syms xVar
syms a

% plateLength
plateLength = 50;

% Nodal coordinates
xNodes = [0, 50, 75, 12.5];
yNodes = [0, 0, 50, 50];
% xNodes = [0, 2, 2, 0];
% yNodes = [0, 0, 2, 2];
% xNodes = [0, a, 3/2*a, 0];
% yNodes = [0, 0, 2, 2*a];


% xNodes = [x1 x2 x3 x4];
% yNodes = [y1 y2 y3 y4];

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

% displacement
w = N1*w1 + N2*w2 + N3*w3 + N4*w4;

wAna = -1/2*MEI*xVar^2+MEI*x1*xVar-phi*xVar+omega+phi*x1-1/2*MEI*x1^2;
% wSubs = simplify(subs(w, [w1 w2 w3 w4], omega + [subs(wAna, xVar, x1) (xNodes(2)^2/plateLength-xNodes(2))*phi (xNodes(3)^2/plateLength-xNodes(3))*phi 0]));
wSubs = simplify(subs(w, [w1 w2 w3 w4], [subs(wAna, xVar, x1) subs(wAna, xVar, x2) subs(wAna, xVar, x3) subs(wAna, xVar, x4)]));
wSubs = simplify(subs(wSubs, [x1 x2 x3 x4], xNodes));
w1Subs = simplify(subs(subs(wAna, xVar, x1), [x1 x2 x3 x4], xNodes));
w2Subs = simplify(subs(subs(wAna, xVar, x2), [x1 x2 x3 x4], xNodes));
w3Subs = simplify(subs(subs(wAna, xVar, x3), [x1 x2 x3 x4], xNodes));
w4Subs = simplify(subs(subs(wAna, xVar, x4), [x1 x2 x3 x4], xNodes));

% Rotations
thetaX = N1*thetaX1 + N2*thetaX2 + N3*thetaX3 + N4*thetaX4;
thetaY = N1*thetaY1 + N2*thetaY2 + N3*thetaY3 + N4*thetaY4;

thetaYAna = -(-MEI*xVar+MEI*x1-phi);
thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
% thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [phi (-2*xNodes(2)/plateLength+1)*phi (-2*xNodes(3)/plateLength+1)*phi phi]));
thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [subs(thetaYAna, xVar, x1) subs(thetaYAna, xVar, x2) subs(thetaYAna, xVar, x3) subs(thetaYAna, xVar, x4)]));
thetaYSubs = simplify(subs(thetaYSubs, [x1 x2 x3 x4], xNodes));

% Jacobian
J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);
J = [J11, J12; J21, J22];

% DSG
int1 = int(diff(wSubs, xi) + J11*thetaYSubs+J21*thetaXSubs, xi);
deltaGamma12 = subs(int1, [xi, eta], [1, -1]) - subs(int1, [xi, eta], [-1, -1]);
deltaGamma13 = subs(int1, [xi, eta], [1, 1]) - subs(int1, [xi, eta], [-1, -1]);
deltaGamma14 = subs(int1, [xi, eta], [-1, 1]) - subs(int1, [xi, eta], [-1, -1]);
GammaXiDSG = diff(N2, xi)*deltaGamma12+diff(N3, xi)*deltaGamma13+diff(N4, xi)*deltaGamma14;
GammaXiDSG = simplify(GammaXiDSG);

int2 = int(diff(wSubs, eta) + J22*thetaXSubs+J12*thetaYSubs, eta);
deltaGamma22 = subs(int2, [xi, eta], [1, -1]) - subs(int2, [xi, eta], [-1, -1]);
deltaGamma23 = subs(int2, [xi, eta], [1, 1]) - subs(int2, [xi, eta], [-1, -1]);
deltaGamma24 = subs(int2, [xi, eta], [-1, 1]) - subs(int2, [xi, eta], [-1, -1]);
GammaEtaDSG = diff(N2, eta)*deltaGamma22+diff(N3, eta)*deltaGamma23+diff(N4, eta)*deltaGamma24;
GammaEtaDSG = simplify(GammaEtaDSG);

% Transverse shear strain
GammaXi = diff(wSubs, xi) + J11*thetaYSubs+J21*thetaXSubs;
GammaEta = diff(wSubs, eta) + J22*thetaXSubs+J12*thetaYSubs;

GammaXi = simplify(collect(simplify(collect(collect(collect(GammaXi, omega), MEI), phi)), [xi, eta]));

% ANS solution
GammaXiANS = 1/2*(1+eta)*subs(GammaXi, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXi, [xi, eta], [0, -1]);
GammaEtaANS = 1/2*(1+xi)*subs(GammaEta, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEta, [xi, eta], [-1, 0]);

% Backtransformation
Gamma = J.'\[GammaXi; GammaEta];
GammaANS = J.'\[GammaXiANS; GammaEtaANS];

disp(collect(Gamma, phi));
disp(GammaANS);

disp(eval(subs(GammaANS, [xi, eta, phi], [-0.5774, -0.5774, 1])));

% Conclusion:
% The ANS interpolation eliminates transverse shear locking for any mesh!
% This is already achieved in the reference element due to the special
% interpolation.

% Recreation:
GammaXiTest = subs(GammaXi, phi, 1);
GammaEtaTest = subs(GammaEta, phi, 1);
% disp('Start solving...');
% sol = solve(GammaXi, [xi, eta]);
% disp(sol);

% There are multiple solutions to this equation!
% Therefore, plotting is necessary
figure;
fplot(subs(GammaXiTest, xi, 0));

sol = solve(subs(GammaXiTest, xi, 0), eta);
disp(sol);


% Final ANS solution
GammaXi = subs(diff(w, xi) + J11*thetaY+J21*thetaX, [xi, eta], [0, 1]);
