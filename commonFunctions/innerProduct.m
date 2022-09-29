function [C] = innerProduct(A, B)
%INNERPRODUCT Creates inner product of two matrices
%
%   Syntax:
%   A, B matrix
%
%   [C] = innerProduct(A,B)
%
% Description
%
% Inner product of matrices.

if size(A) == size(B)
    C = trace(A.'*B);
else
    error('Matrix dimensions do not agree!')
end
end

