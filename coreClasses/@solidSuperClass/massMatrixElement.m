function dataFE = massMatrixElement(obj,setupObject)
% compute constant mass matrix
%% setup
% load objects
meshObject = obj.meshObject;
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
%   mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof,1);
dimension = obj.dimension;
additionalFields = obj.additionalFields;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% nodal dofs
qR = meshObject.nodes;
% material data and voigt notation
rho = materialObject.rho;
% initialization energy and storageFEObject
[dataFE, initialDataFE] = storageFEObject.initializeMassDataFE;
%% element loop for residual and tangent
parfor e = 1:numberOfElements
    numberOfDofs = numel(globalFullEdof(e,:));
    Me = zeros(numberOfDofs);
    edR = qR(edof(e, :), 1:dimension)';    
    for k = 1:numberOfGausspoints
        dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
        % compute Jacobi matrix and determinant
        [~,detJ] = computeJacobian(edR,dN_xi_I,setupObject.toleranceDetJ);
%         Nkron = kron(N(k,:)',I);
%         index2 = 1:numberOfDofs;
%         if additionalFields > 0
%             index2(dimension+additionalFields:dimension+additionalFields:numberOfDofs) = [];
%         end
%         Me(index2,index2)  = Me(index2,index2) + rho*(Nkron*Nkron')*detJ*gaussWeight(k);
        A1 = (N_k_I(k,:)'*N_k_I(k,:))*rho;
        MAT = zeros(numberOfDofs);
        sizeN = size(N_k_I,2)*(dimension + additionalFields);
        for l = 1:dimension
            MAT(l:dimension+additionalFields:sizeN,l:dimension+additionalFields:sizeN) = A1;
        end
        Me = Me + MAT*detJ*gaussWeight(k);
    end
    elementDataFE = storageFEObject.assignElementMassDataFE(initialDataFE,Me(:),globalFullEdof,e);
    dataFE(e) = elementDataFE;
end
end
