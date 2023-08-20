function [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension)
%EXTRACTJACOBIANFORGAUSSPOINT Summary of this function goes here
%   Detailed explanation goes here

J = JAll(:, dimension*(k - 1)+1:dimension*k);
if size(J, 1) == size(J, 2)
    detJ = det(J);
elseif size(J, 2) == 1 && size(J, 1) ~= 1
    detJ = norm(J);
else
    error('not implemented yet.')
end

if setupObject.computePostData
    detJ = abs(detJ); % Jacobian determinant
end

% check the Jacobian determinant
if detJ < 10 * setupObject.toleranceDetJ
    error('Jacobi determinant equal or less than zero.')
end
end
