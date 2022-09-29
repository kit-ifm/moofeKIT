function F0 = F0Matrix(dimension, J0)
% F0MATRIX Returns F0 matrix
% This function returns the F0 matrix, which is used in connection with
% Pian-Sumihara and EAS formulations.
% For Pian-Sumihara formulations the F0 matrix results from the
% transformation of the ansatz for the stress
% sigma = J0 * Sigma * J0'
% to it's equivalent in Voigt's vector notation
% sigmaVoigt = F0 * SigmaVoigt.
% For EAS formulations the F0 matrix results from the transformation of the
% ansatz for the enhanced part of the strain field
% epsilonTilde = j0/j * inv(J0)' * ETilde * inv(J0)
% to it's equivalent in Voigt's vector notation
% epsilonTildeVoigt = j0/j * inv(F0)' * ETildeVoigt.
% The F0 matrix is the same for both Pian-Sumihara and EAS formulations.
%
% CALL
% F0 = F0Matrix(dimension, J0)
% dimension: number of dimensions
% J0: Jacobian matrix evaluated at the center of the element
% F0: F0 matrix
%
% REFERENCE
% Lecture notes FEidF
%
% CREATOR(S)
% Felix Zaehringer

if dimension == 2
    F0 = [J0(1, 1) * J0(1, 1), J0(1, 2) * J0(1, 2), 2 * J0(1, 1) * J0(1, 2); ...
        J0(2, 1) * J0(2, 1), J0(2, 2) * J0(2, 2), 2 * J0(2, 1) * J0(2, 2); ...
        J0(1, 1) * J0(2, 1), J0(1, 2) * J0(2, 2), J0(1, 1) * J0(2, 2) + J0(1, 2) * J0(2, 1)];
elseif dimension == 3
    F0 = [J0(1, 1)^2, J0(1, 2)^2, J0(1, 3)^2, 2 * J0(1, 1) * J0(1, 2), 2 * J0(1, 2) * J0(1, 3), 2 * J0(1, 1) * J0(1, 3); ...
        J0(2, 1)^2, J0(2, 2)^2, J0(2, 3)^2, 2 * J0(2, 1) * J0(2, 2), 2 * J0(2, 2) * J0(2, 3), 2 * J0(2, 1) * J0(2, 3); ...
        J0(3, 1)^2, J0(3, 2)^2, J0(3, 3)^2, 2 * J0(3, 1) * J0(3, 2), 2 * J0(3, 2) * J0(3, 3), 2 * J0(3, 1) * J0(3, 3); ...
        J0(1, 1) * J0(2, 1), J0(1, 2) * J0(2, 2), J0(1, 3) * J0(2, 3), J0(1, 1) * J0(2, 2) + J0(1, 2) * J0(2, 1), J0(1, 2) * J0(2, 3) + J0(1, 3) * J0(2, 2), J0(1, 1) * J0(2, 3) + J0(1, 3) * J0(2, 1); ...
        J0(2, 1) * J0(3, 1), J0(2, 2) * J0(3, 2), J0(2, 3) * J0(3, 3), J0(2, 1) * J0(3, 2) + J0(2, 2) * J0(3, 1), J0(2, 2) * J0(3, 3) + J0(2, 3) * J0(3, 2), J0(2, 1) * J0(3, 3) + J0(2, 3) * J0(3, 1); ...
        J0(1, 1) * J0(3, 1), J0(1, 2) * J0(3, 2), J0(1, 3) * J0(3, 3), J0(1, 1) * J0(3, 2) + J0(1, 2) * J0(3, 1), J0(1, 2) * J0(3, 3) + J0(1, 3) * J0(3, 2), J0(1, 1) * J0(3, 3) + J0(1, 3) * J0(3, 1)];
else
    error('F0 matrix is not implemented for this dimension');
end

end
