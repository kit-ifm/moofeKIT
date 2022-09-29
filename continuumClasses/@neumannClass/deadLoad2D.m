function [rData, kData, elementEnergy, array] = deadLoad2D(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

% load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dimension = obj.masterObject.dimension - 1;
I = eye(dimension+1);
edofNo = meshObject.edof;
qR = obj.masterObject.qR;

switch lower(obj.projectionType)
    case 'none'
        projection = 0;
    case 'x'
        projection = 1;
    case 'y'
        projection = 2;
    otherwise
        error('type of projection not implemented')
end
forceVector = obj.forceVector * obj.timeFunction(obj.time);

elementEnergy.externalEnergy = 0;
edofLocal = edofNo(e, :);
edR = qR(edofLocal, 1:dimension+1)';
for k = 1:numberOfGausspoints
    N_I = N(k, :);
    dN_xiI = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
    switch projection
        case 0
            J = edR * dN_xiI';
        case 1
            J = edR(1, :) * dN_xiI';
        case 2
            J = edR(2, :) * dN_xiI';
    end
    detJ = norm(J);
    if detJ <= 10 * eps
        error('Jacobi-Determinant less or equal to zero')
    end
    rData{1} = rData{1} - kron(N_I, I)' * forceVector * detJ * gaussWeight(k);
end

end