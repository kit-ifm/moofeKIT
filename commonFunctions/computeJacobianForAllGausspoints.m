function JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I)
%COMPUTEJACOBIANFORALLGAUSSPOINTS Summary of this function goes here
%   Detailed explanation goes here

dNr = reshape(dN_xi_k_I, [], size(dN_xi_k_I, 3));
JAll = ed * dNr';
end
