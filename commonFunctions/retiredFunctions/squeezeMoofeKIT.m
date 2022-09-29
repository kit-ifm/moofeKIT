function B = squeezeMoofeKIT(A)
%SQUEEZE Remove singleton dimensions.
%   B = SQUEEZE(A) returns an array B with the same elements as
%   A but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(A,dim)==1.  2-D arrays are
%   unaffected by squeeze so that row vectors remain rows.
%
%   For example,
%       squeeze(rand(2,1,3))
%   is 2-by-3.
%
%   See also SHIFTDIM.
%   Copyright 1984-2020 The MathWorks, Inc.
%   Edited by Marlon Franke
%   compare with B = reshape(A(:,k,:),[size(A,1),size(A,3)])

if ~ismatrix(A)
    sizeA = size(A);       % siz has at least one element not equal to 1
    sizeA(sizeA == 1) = [];  % remove singleton dimensions
    if isscalar(sizeA)
        B = reshape(A,1,sizeA);  % reshape to 2-D
    else
        B = reshape(A,sizeA);
    end
else
    B = A;
end
