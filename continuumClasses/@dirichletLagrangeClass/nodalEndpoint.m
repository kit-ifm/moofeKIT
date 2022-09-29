function nodalEndpoint(obj,setupObject,varargin)
% Elementroutine nodwise dirichlet with lagrange multiplicator
globalFullEdof = obj.meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof,1);
numberOfDofs = size(globalFullEdof,2);
qR = obj.masterObject.qR(obj.nodeList,obj.nodalDof);
qN1 = obj.masterObject.qN1(obj.nodeList,obj.nodalDof);
lambdaN1 = obj.qN1;
%prescribed displacements
uDirichlet = obj.timeFunction(obj.time);
storageFEObject = obj.storageFEObject;
initializeDataFE(storageFEObject);
for e = 1:numberOfElements
    Re = zeros(numberOfDofs,1);
    Ke = zeros(numberOfDofs,numberOfDofs);
    % dofs
    q = qN1(e);
    qDirichelt = qR(e) + uDirichlet;
    lambda = lambdaN1(e);
    % residual
    Re(1) = -lambda;
    Re(2) = -(q-qDirichelt);
    % tangent
    Ke(1,2) = -1;
    Ke(2,1) = -1;    
    % output data
    storeDataFE(storageFEObject,Re,Ke(:),globalFullEdof,e);
end

