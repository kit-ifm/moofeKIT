%% Quadrilateral with J0
syms xi eta;
syms x1 x2 x3 x4;
syms y1 y2 y3 y4;
syms w1 w2 w3 w4;
syms thetaX1 thetaX2 thetaX3 thetaX4;
syms thetaY1 thetaY2 thetaY3 thetaY4;
syms alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7 alpha8;

N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);

x = N1*x1+N2*x2+N3*x3+N4*x4;
y = N1*y1+N2*y2+N3*y3+N4*y4;


J = [diff(x, xi), diff(x, eta); diff(y, xi), diff(y, eta)];
% J = subs(J, [xi, eta], [0, 0]);
J11 = J(1, 1);
J12 = J(1, 2);
J21 = J(2, 1);
J22 = J(2, 2);
detJ = det(J);

J0 = subs(J, {xi, eta}, {0, 0});
J011 = J0(1, 1);
J012 = J0(1, 2);
J021 = J0(2, 1);
J022 = J0(2, 2);
detJ0 = det(J0);

JInv = inv(J.');
J0Inv = inv(J0.');

w = N1*w1+N2*w2+N3*w3+N4*w4;
thetaX = N1*thetaX1+N2*thetaX2+N3*thetaX3+N4*thetaX4;
thetaY = N1*thetaY1+N2*thetaY2+N3*thetaY3+N4*thetaY4;

GammaXiO = diff(w, xi)+J011*thetaX+J021*thetaY;
GammaEtaO = diff(w, eta)+J012*thetaX+J022*thetaY;
[GammaXiC, GammaXiT] = coeffs(GammaXiO, [xi, eta]);
[GammaEtaC, GammaEtaT] = coeffs(GammaEtaO, [xi, eta]);
GammaXi = GammaXiC*(GammaXiT.*[alpha1, alpha2, alpha3, alpha4]).';
GammaEta = GammaEtaC*(GammaEtaT.*[alpha5, alpha6, alpha7, alpha8]).';

GammaX = simplify(JInv(1, :) * [GammaXi; GammaEta]);
GammaY = simplify(JInv(2, :) * [GammaXi; GammaEta]);

% Bending
GammaXBending = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, 1/2*x4^2+x4+1-y4^2, -x1-1, -x2-1, -x3-1, -x4-1, 2*y1, 2*y2, 2*y3, 2*y4]);
GammaYBending = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, 1/2*x4^2+x4+1-y4^2, -x1-1, -x2-1, -x3-1, -x4-1, 2*y1, 2*y2, 2*y3, 2*y4]); %[1/2*x1^2+x1+1-y1^2, 1/2*x2^2+x2+1-y2^2, 1/2*x3^2+x3+1-y3^2, -x1-1, -x2-1, -x3-1, 2*y1, 2*y2, 2*y3]

GammaXBending2 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);
GammaYBending2 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1/2*y1^2, 1/2*y2^2, 1/2*y3^2, 1/2*y4^2, 0, 0, 0, 0, -y1, -y2, -y3, -y4]);

% Shearing
GammaXShearing = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1, x2, x3, x4, 0, 0, 0, 0, 0, 0, 0, 0]);

GammaXShearing2 = subs(GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);
GammaYShearing2 = subs(GammaY, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1, y2, y3, y4, 0, 0, 0, 0, 0, 0, 0, 0]);

% solving
sol = solve(subs([GammaXBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaYBending; GammaYBending2; GammaYShearing; GammaYShearing2], [x3, x4, y2, y3], [x2, x1, y1, y4]) == [0; 0; 1; 0; 0; 0; 0; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
solO = solve([GammaXBending; GammaXBending2; GammaXShearing; GammaXShearing2; GammaYBending; GammaYBending2; GammaYShearing; GammaYShearing2] == [0; 0; 1; 0; 0; 0; 0; 1], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]);
GammaXiCorrect = simplify(subs(GammaXi, [alpha1, alpha2, alpha3, alpha4], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4]));
GammaEtaCorrect = simplify(subs(GammaEta, [alpha5, alpha6, alpha7, alpha8], [sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]));

GammaXCorrect = simplify(subs(GammaX, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]));
GammaYCorrect = simplify(subs(GammaY, [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8], [sol.alpha1, sol.alpha2, sol.alpha3, sol.alpha4, sol.alpha5, sol.alpha6, sol.alpha7, sol.alpha8]));

ANSSolution_GammaXi = 1/2*(1+eta)*subs(GammaXiO, [xi, eta], [0, 1])+1/2*(1-eta)*subs(GammaXiO, [xi, eta], [0, -1]);
ANSSolution_GammaEta = 1/2*(1+xi)*subs(GammaEtaO, [xi, eta], [1, 0])+1/2*(1-xi)*subs(GammaEtaO, [xi, eta], [-1, 0]);
ANSSolution_GammaX = J0Inv(1, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];
ANSSolution_GammaY = J0Inv(2, :) * [ANSSolution_GammaXi; ANSSolution_GammaEta];

syms absx absy;
subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [1,1,1,1,0,0,0,0,0,0,0,0]);
subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1,x2,x3,x4,0,0,0,0,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]);
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1,y2,y3,y4,0,0,0,0,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX-y, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1,x2*y2,x3*y3,x4*y4,0,0,0,0,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]))
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,1,1,1,1,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,0,0,0,0,1,1,1,1]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,x1,x2,x3,x4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,0,0,0,0,x1,x2,x3,x4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX-y, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,y1,y2,y3,y4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,0,0,0,0,y1,y2,y3,y4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX-x*y, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,x1*y1,x2*y2,x3*y3,x4*y4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,0,0,0,0,x1*y1,x2*y2,x3*y3,x4*y4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaXi, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1,x2*y2,x3*y3,x4*y4,-y1,-y2,-y3,-y4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1,x2*y2,x3*y3,x4*y4,-y1,-y2,-y3,-y4,-x1,-x2,-x3,-x4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2,x2^2,x3^2,x4^2,-2*x1,-2*x2,-2*x3,-2*x4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [y1^2,y2^2,y3^2,y4^2,0,0,0,0,-2*y1,-2*y2,-2*y3,-2*y4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2*y1,x2^2*y2,x3^2*y3,x4^2*y4,-2*x1*y1,-2*x2*y2,-2*x3*y3,-2*x4*y4,-x1^2,-x2^2,-x3^2,-x4^2]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^3,x2^3,x3^3,x4^3,-3*x1^2,-3*x2^2,-3*x3^2,-3*x4^2,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2,x2^2,x3^2,x4^2,-2*x1,-2*x2,-2*x3,-2*x4,x1,x2,x3,x4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1,x2,x3,x4,-1,-1,-1,-1,1,1,1,1]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2,x2^2,x3^2,x4^2,-2*x1,-2*x2,-2*x3,-2*x4,1,1,1,1]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2*y1,x2^2*y2,x3^2*y3,x4^2*y4,-2*x1*y1,-2*x2*y2,-2*x3*y3,-2*x4*y4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX-2*x, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2,x2^2,x3^2,x4^2,0,0,0,0,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaX, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1^2*y1^2,x2^2*y2^2,x3^2*y3^2,x4^2*y4^2,-2*x1*y1^2,-2*x2*y2^2,-2*x3*y3^2,-2*x4*y4^2,-2*x1^2*y1,-2*x2^2*y2,-2*x3^2*y3,-2*x4^2*y4]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
simplify(subs(subs(ANSSolution_GammaXi-J011*x, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [0,0,0,0,x1,x2,x3,x4,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14]));
subs(subs(ANSSolution_GammaX-y, [w1, w2, w3, w4, thetaX1, thetaX2, thetaX3, thetaX4, thetaY1, thetaY2, thetaY3, thetaY4], [x1*y1,x2*y2,x3*y3,x4*y4,0,0,0,0,0,0,0,0]), [x1, x2, x3, x4, y1, y2, y3, y4], [1, 8, 9, 2, 1, 3, 16, 14])


differenceGammaX = simplify(ANSSolution_GammaX - GammaXCorrect)
differenceGammaY = simplify(ANSSolution_GammaY - GammaYCorrect)

