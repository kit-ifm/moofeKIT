function [rData, kData, elementEnergy, array] = deadLoadMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% The reference volume is used, bodyforce may vary along position
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

if isa(obj.masterObject, 'stringClass') || isa(obj.masterObject, 'beamClass')
    dimension = size(obj.masterObject.meshObject.nodes,2);
else
    dimension = obj.masterObject.dimension;
end
edof = obj.masterObject.meshObject.edof;
DT = setupObject.timeStepSize;
timeFunction = obj.timeFunction(obj.time- 0.5*DT);
loadFunction = obj.loadFunction;

if isa(obj.masterObject, 'beamClass')
    edR = obj.masterObject.meshObject.nodes(edof(e,:), 1:dimension)';
else
    edR = obj.masterObject.qR(edof(e,:), 1:dimension)';
end

edN = obj.masterObject.qN(edof(e,:), 1:dimension)';
edN1 = obj.masterObject.qN1(edof(e,:), 1:dimension)';
edN05 = 1/2*(edN+ edN1);

elementEnergy.externalEnergy = 0;
for k = 1:numberOfGausspoints
    dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
    [~,detJ] = computeJacobian(edR,dN_xi_I,setupObject.toleranceDetJ,setupObject.computePostData);
    %position and bodyforce
    Xh = kron(N_k_I(k,:),eye(dimension))*edR(:);
    xh = kron(N_k_I(k,:),eye(dimension))*edN1(:);
    if isa(loadFunction,'function_handle')
        B0 = loadFunction(Xh(1),Xh(2),Xh(3));
    else
        B0 = loadFunction;
    end
    %residual
    if isa(obj.masterObject, 'beamClass')
        rData{1} = rData{1} - timeFunction*kron(N_k_I(k,:),eye(obj.dimension))'*B0*detJ*gaussWeight(k);
    else
        rData{1} = rData{1} - timeFunction*kron(N_k_I(k,:),eye(dimension))'*B0*detJ*gaussWeight(k);
    end
    
    elementEnergy.externalEnergy = elementEnergy.externalEnergy - xh'*B0*detJ*gaussWeight(k);
end
end