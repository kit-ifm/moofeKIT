% Script to investigate transverse shear locking for Bubnov-Galerkin FE

% Syms
syms xi eta
syms x1 y1 x2 y2 x3 y3 x4 y4
syms thetaX1 thetaX2 thetaX3 thetaX4
syms thetaY1 thetaY2 thetaY3 thetaY4
syms phi

% bendingType: bending around which axis
bendingType = 'y';

% elementForm
distortionType = 'angular';

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

% Rotations
thetaX = N1*thetaX1 + N2*thetaX2 + N3*thetaX3 + N4*thetaX4;
thetaY = N1*thetaY1 + N2*thetaY2 + N3*thetaY3 + N4*thetaY4;

if strcmp(bendingType, 'x')
    thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [phi phi -phi -phi]));
    thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [0 0 0 0]));
elseif strcmp(bendingType, 'y')
    thetaXSubs = simplify(subs(thetaX, [thetaX1 thetaX2 thetaX3 thetaX4], [0 0 0 0]));
    thetaYSubs = simplify(subs(thetaY, [thetaY1 thetaY2 thetaY3 thetaY4], [phi -phi -phi phi]));
end

% Jacobian
J11 = diff(x, xi);
J12 = diff(x, eta);
J21 = diff(y, xi);
J22 = diff(y, eta);

% Transverse shear strain
GammaX = J11*thetaYSubs+J21*thetaXSubs;
GammaY = J22*thetaXSubs+J12*thetaYSubs;

% Find matching xi/eta for Kirchhoff hypothesis


disp('Start solving...');
sol = solve([GammaX, GammaY], [xi, eta]);
disp(GammaX);
disp(GammaY);
disp(sol);

% Conclusion:
% The Kirchhoff condition is satisfied for any element geometry, when
% interpolating GammaX at xi=0 (eta is arbitrary) and GammaY at eta=0 (xi
% is arbitrary).

% Correct solution
GammaX = subs(J11*thetaY+J21*thetaX, xi, 0);
GammaY = subs(J22*thetaX+J12*thetaY, eta, 0);
disp(GammaX);
disp(GammaY);

% ANS solution
GammaXANS = 1/2*(1+eta)*subs(J11*thetaY+J21*thetaX, [xi, eta], [0, 1])+1/2*(1-eta)*subs(J11*thetaY+J21*thetaX, [xi, eta], [0, -1]);

disp(GammaXANS);
disp(simplify(GammaX-GammaXANS))

