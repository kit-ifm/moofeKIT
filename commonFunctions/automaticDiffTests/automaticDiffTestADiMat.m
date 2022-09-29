ADiMat_startup

%% 1
% A = magic(2);
% c = [ 2 3 4 5];
% f2(A, c)

% Jac = admDiffFor(@f, 1, A, c);
% Jac = admDiffFor(@f, 1, A, c, admOptions('flags', '--check-certificate'));
% Jac = admDiffRev(@f, 1, A, c, admOptions('flags', '--check-certificate'));

%% 2
% x = [1 2 3 4];
% y = rosenbrock(x)
% Jac = admDiffRev(@rosenbrock, 1, x, admOptions('flags', '--check-certificate'));

%% 3
x = [1 2 3 4];
a = 3; b = 4; c = 5;
y = simpleFunction(x,a,b,c)
[Jac y2] = admDiffRev(@simpleFunction, 1, x, a, b,c, admOptions('flags', '--check-certificate'))
dFunction = 2*x.*a.*b
[Jac2 b2 b3 b4] = a_simpleFunction(x,a,b,c,eye(4))

%% 4
% x = [1 2 3 4];
% a = 3; b = 4; c = 5;
% y = simpleFunction(x,a,b,c)
% [Jac y2] = admDiffFor(@simpleFunction, 1, x, a, b,c, admOptions('flags', '--check-certificate'))
% dFunction = 2*x.*a.*b
% [g_y, y]= g_simpleFunction(1, x, 0, a, 0, b, 0, c)
