function [J,detJ] = computeJacobian(ed, dN_xi_I, tolerance, varargin)
% COMPUTEJACOBIAN Function to compute the Jacobian determinant
%
% CALL [J, detJk] = computeJacobian(ed,dNxik,tolerance)
% input:    ed denotes the geometric nodes of the considered element
%           dN_xi_I denotes the derivative of the shapefunctions wrt to xi at k-th gausspoint
%           tolerance of the Jacobian determinant
% output:   J denotes the Jacobian matrix
%           detJ Jacobian determinant for reference configuration

J = ed * dN_xi_I';          % Jacobian matrix
if size(J, 1) == size(J, 2)
    detJ = det(J);
elseif size(J, 2) == 1 && size(J, 1) ~= 1
    detJ = norm(J);
else
    error('not implemented yet.')
end

if numel(varargin) >= 1
    computePostData = varargin{1};
    if computePostData
        detJ = abs(det(J));              % Jacobian determinant
    end
end

% check the Jacobian determinant
if detJ < 10 * tolerance
    error('Jacobi determinant equal or less than zero.')
end
end