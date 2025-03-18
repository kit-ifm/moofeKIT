%% Triangle
syms xi eta;
syms x1 x2 x3;
syms y1 y2 y3;
syms w1 w2 w3;
syms thetaX1 thetaX2 thetaX3;
syms thetaY1 thetaY2 thetaY3;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6;

N1 = 1-xi-eta;
N2 = xi;
N3 = eta;

x = N1*x1+N2*x2+N3*x3;
y = N1*y1+N2*y2+N3*y3;

J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);
J = [J11, J12; J21, J22];
detJ = det(J);

JInv = inv(J.');

w = N1*w1+N2*w2+N3*w3;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3;

% GammaXi = coeffs(diff(w, xi)+J11*thetaX+J21*thetaY, [xi, eta])*[alpha1; alpha2*xi; alpha3*eta];
% GammaEta = coeffs(diff(w, eta)+J12*thetaX+J22*thetaY, [xi, eta])*[alpha4; alpha5*xi; alpha6*eta];

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;
[GammaXiC, GammaXiT] = coeffs(GammaXiO, [xi, eta]);
[GammaEtaC, GammaEtaT] = coeffs(GammaEtaO, [xi, eta]);
GammaXi = GammaXiC*(GammaXiT.*[alpha1, alpha2, alpha3]).';
GammaEta = GammaEtaC*(GammaEtaT.*[alpha4, alpha5, alpha6]).';

GammaX = simplify(JInv(1, :) * [GammaXi; GammaEta]);
GammaY = simplify(JInv(2, :) * [GammaXi; GammaEta]);

% Bending
GammaXBending = subs(GammaX, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, -x1-1, -x2-1, -x3-1, 2*y1, 2*y2, 2*y3]);
GammaYBending = subs(GammaY, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, -x1-1, -x2-1, -x3-1, 2*y1, 2*y2, 2*y3]); %[1/2*x1^2, 1/2*x2^2, 1/2*x3^2, -x1, -x2, -x3, 0, 0, 0]

GammaXBending2 = subs(GammaX, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 0, 0, 0, -y1, -y2, -y3]);
GammaYBending2 = subs(GammaY, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 0, 0, 0, -y1, -y2, -y3]);

% Shearing
GammaXShearing = subs(GammaX, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [x1, x2, x3, 0, 0, 0, 0, 0, 0]);
GammaYShearing = subs(GammaY, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [x1, x2, x3, 0, 0, 0, 0, 0, 0]);

GammaXShearing2 = subs(GammaX, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [y1, y2, y3, 0, 0, 0, 0, 0, 0]);
GammaYShearing2 = subs(GammaY, [w1, w2, w3, thetaX1, thetaX2, thetaX3, thetaY1, thetaY2, thetaY3], [y1, y2, y3, 0, 0, 0, 0, 0, 0]);

% solving
sol = solve([GammaXBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaYBending; GammaYBending2; GammaYShearing; GammaYShearing2] == [0; 0; 1; 0; 0; 0; 0; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6]);
GammaXiCorrect = simplify(subs(GammaXi, [alpha1, alpha2, alpha3], [sol.alpha1, sol.alpha2, sol.alpha3]));
GammaEtaCorrect = simplify(subs(GammaEta, [alpha4, alpha5, alpha6], [sol.alpha4, sol.alpha5, sol.alpha6]));

GammaXCorrect = simplify(subs(GammaX, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6]));
GammaYCorrect = simplify(subs(GammaY, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6]));

DSGSolution_GammaX = JInv(1, :) * [(w2-w1+1/2*J11*(thetaX1 + thetaX2)+1/2*J21*(thetaY1 + thetaY2)); (w3-w1+1/2*J12*(thetaX1 + thetaX3)+1/2*J22*(thetaY1 + thetaY3))];
DSGSolution_GammaY = JInv(2, :) * [(w2-w1+1/2*J11*(thetaX1 + thetaX2)+1/2*J21*(thetaY1 + thetaY2)); (w3-w1+1/2*J12*(thetaX1 + thetaX3)+1/2*J22*(thetaY1 + thetaY3))];

differenceGammaX = simplify(DSGSolution_GammaX - GammaXCorrect)
differenceGammaY = simplify(DSGSolution_GammaY - GammaYCorrect)

