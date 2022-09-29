function [x] = solverLinearSystem(setupObject, A, b)
%SOLVERLINEARSYSTEM Returns the solution to a system of linear equations
%
% CALL: x = solverLinearSystem(setupObject, A, b)
% setupObject: setupObject
% A: coefficient matrix
% b: right side vector
%
% OUTPUT:
% x: solution vector

% check input
assert(size(A, 1) == size(A, 2), 'Coefficient matrix must be square!');
assert(size(A, 1) == size(b, 1), 'Coefficient matrix and right side must have the same number of rows!');
assert(~any(any(isnan(b))), "Right side vector contains invalid entries!");
assert(~any(any(isnan(A))), "Coefficient matrix contains invalid entries!");


% check if preconditioning should be done
flagPreconditioning = false;
if setupObject.usePreconditioning
    if cond(A) > 1e10
        flagPreconditioning = true;
        [P, R, C] = equilibrate(A);
        A = R * P * A * C;
        b = R * P * b;
    end
end


% solve linear system
x = A \ b;

if flagPreconditioning
    x = C * x;
end
end
