function [rData, kData, elementEnergy, array] = deadLoadEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% The reference volume is used, bodyforce may vary along position
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dimension = obj.masterObject.dimension;
edof = obj.masterObject.meshObject.edof;
timeFunction = obj.timeFunction(obj.time);
loadFunction = obj.loadFunction;
edR = obj.masterObject.qR(edof(e,:), 1:dimension)';
elementEnergy.externalEnergy = 0;
for k = 1:numberOfGausspoints
    dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
    [~,detJ] = computeJacobian(edR,dN_xi_I,setupObject.toleranceDetJ,setupObject.computePostData);
    %position and bodyforce
    Xh = kron(N_k_I(k,:),eye(dimension))*edR(:);
    if isa(loadFunction,'function_handle')
        B0 = loadFunction(Xh(1),Xh(2),Xh(3));
    else
        B0 = loadFunction;
    end
    %residual
    rData{1} = rData{1} - timeFunction*kron(N_k_I(k,:),eye(dimension))'*B0*detJ*gaussWeight(k);
    elementEnergy.externalEnergy = elementEnergy.externalEnergy - xh'*B0*detJ*gaussWeight(k);
end
end