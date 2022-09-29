function [rData, kData, elementEnergy, array] = deadLoad(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
masterObject = obj.masterObject;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

dimension = obj.masterObject.dimension - 1;
edof = obj.meshObject.edof;
qR = masterObject.qR;
timeStep = setupObject.timeStep;
elementEnergy.externalEnergy = 0;

if strcmpi(obj.field,'mechanical')
        forceVector = obj.forceVector * obj.timeFunction(obj.time);
        edofLocal = edof(e, :);
        edR = qR(edofLocal, 1:dimension+1)';
        for k = 1:numberOfGausspoints
%             XREF = kron(N_k_I(k, :), eye(dimension+1)) * edR(:);
            dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
            J = edR * dN_xi_I';
            detA0 = norm(cross(J(:,1),J(:,2)));
            rData{1} = rData{1} - kron(N_k_I(k, :)', forceVector) * detA0 * gaussWeight(k);
        end
elseif strcmpi(obj.field,'electrical')
        QN1 = masterObject.qN1(:,1:obj.masterObject.dimension+1);
        QREF2 = masterObject.qR(:,1:obj.masterObject.dimension+1);
        omega_0 = obj.omega_0*obj.timeFunction(obj.time);
        edofLocal = edof(e,:);
        phiN1 = QN1(edofLocal,dimension+1)';
        phiRef = QREF2(edofLocal,dimension+1)';
        ePot = 0;
        edR = qR(edofLocal, 1:dimension+1)';
        for k = 1:numberOfGausspoints
            dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
            J = edR * dN_xi_I';
            detA0 = norm(cross(J(:,1),J(:,2)));
            rData{1} = rData{1} + N_k_I(k,:)'*omega_0*detA0* gaussWeight(k);
            ePot = ePot + ((phiN1-phiRef)*N_k_I(k,:)')*omega_0*detA0* gaussWeight(k);
        end
elseif strcmpi(obj.field,'thermal')
        energy = obj.energy*obj.timeFunction(obj.time);
        edR = qR(edof(e,:),1:dimension+1)';
        for k = 1:numberOfGausspoints
            dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
            J = edR * dN_xi_I';
            detA0 = norm(cross(J(:,1),J(:,2)));
%             detA0 = norm(cross(JR(:,1+(k-1)*2),JR(:,2+(k-1)*2)));
            rData{1} = rData{1} - N_k_I(k,:)'*energy*detA0*gaussWeight(k);
        end
else
    error('Field not implemented')
end
end