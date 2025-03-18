%% Quadrilateral
syms xi eta;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8 alpha9 alpha10 alpha11 alpha12;

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

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J11*thetaX+J21*thetaY;
GammaEtaO = diff(w, eta)+J12*thetaX+J22*thetaY;
[GammaXiC, GammaXiT] = coeffs(GammaXiO, [xi, eta]);
[GammaEtaC, GammaEtaT] = coeffs(GammaEtaO, [xi, eta]);
GammaXi = GammaXiC*(GammaXiT.*[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6]).';
GammaEta = GammaEtaC*(GammaEtaT.*[alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]).';

GammaX = simplify(JInv(1, :) * [GammaXi; GammaEta]);
GammaY = simplify(JInv(2, :) * [GammaXi; GammaEta]);

% Bending
GammaXiBending = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, 1/2*x4^2+x4+1-y4^2, -x1-1, -x2-1, -x3-1, -x4-1, 2*y1, 2*y2, 2*y3, 2*y4]);
GammaEtaBending = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, 1/2*x4^2+x4+1-y4^2, -x1-1, -x2-1, -x3-1, -x4-1, 2*y1, 2*y2, 2*y3, 2*y4]); %[1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, -x1-1, -x2-1, -x3-1, 2*y1, 2*y2, 2*y3]

GammaXiBendingRegular = subs(GammaXiBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBendingRegular = subs(GammaEtaBending, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBendingLinear = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);
GammaEtaBendingLinear = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2*y1, 1/2*x2^2*y2, 1/2*x3^2*y3, 1/2*x4^2*y4, -x1*y1, -x2*y2, -x3*y3, -x4*y4, -1/2*x1^2, -1/2*x2^2, -1/2*x3^2, -1/2*x4^2]);

GammaXiBending2 = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
GammaEtaBending2 = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);

GammaXiBending2Regular = subs(GammaXiBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBending2Regular = subs(GammaEtaBending2, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

GammaXiBendingQuadratic = subs(GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2*y1^2, x2^2*y2^2, x3^2*y3^2, x4^2*y4^2, -2*x1*y1^2, -2*x2*y2^2, -2*x3*y3^2, -2*x4*y4^2, -2*x1^2*y1, -2*x2^2*y2, -2*x3^2*y3, -2*x4^2*y3]);
GammaEtaBendingQuadratic = subs(GammaEta, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2*y1^2, x2^2*y2^2, x3^2*y3^2, x4^2*y4^2, -2*x1*y1^2, -2*x2*y2^2, -2*x3*y3^2, -2*x4*y4^2, -2*x1^2*y1, -2*x2^2*y2, -2*x3^2*y3, -2*x4^2*y3]);

GammaXiBendingQuadraticRegular = subs(GammaXiBendingQuadratic, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);
GammaEtaBendingQuadraticRegular = subs(GammaEtaBendingQuadratic, [x1,x2,x3,x4,y1,y2,y3,y4], [-1,1,1,-1,-1,-1,1,1]);

% Shearing
GammaXShearing = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXShearing2 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing2 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);

% solving
% solA = solve([GammaXiBending; GammaXiBending2; GammaXShearing; GammaXShearing2; GammaEtaBending; GammaEtaBending2; GammaYShearing; GammaYShearing2; alpha1; alpha2; alpha3;alpha7;alpha9;alpha11] == [0; 0; 1; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);
sol = solve([GammaXiBending; GammaXiBending2; GammaXShearing; GammaXShearing2; GammaEtaBending; GammaEtaBending2; GammaYShearing; GammaYShearing2; GammaXiBendingRegular; GammaEtaBendingRegular; GammaXiBending2Regular; GammaEtaBending2Regular] == [0; 0; 1; 0; 0; 0; 0; 1; 0; 0; 0; 0], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12]);
% GammaXiCorrect = simplify(subs(GammaXi, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6]));
% GammaEtaCorrect = simplify(subs(GammaEta, [alpha7, alpha8, alpha9, alpha10, alpha11, alpha12], [sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10, sol.alpha11, sol.alpha12]));
% 
% GammaXCorrect = simplify(subs(GammaX, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10, sol.alpha11, sol.alpha12]));
% GammaYCorrect = simplify(subs(GammaY, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10, alpha11, alpha12], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8, sol.alpha9, sol.alpha10, sol.alpha11, sol.alpha12]));

ANSSolution_GammaXi = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
ANSSolution_GammaEta = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);
ANSSolution_GammaX = JInv(1, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];
ANSSolution_GammaY = JInv(2, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];

differenceGammaX = simplify(ANSSolution_GammaX - GammaXCorrect)
differenceGammaY = simplify(ANSSolution_GammaY - GammaYCorrect)

