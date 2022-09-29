function displacementSCSaintVenantAutomaticDiffEndpoint(obj,setupObject)
%% setup
% load objects
meshObject = obj.meshObject;
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof,1);
numberOfDofs = size(globalFullEdof,2);
dimension = obj.dimension;
% gauss integration and shape function
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
%   N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
qN1 = obj.qN1;
% material data and voigt notation
lambda = materialObject.lambda;
mu = materialObject.mu;
DMat = [lambda+2*mu lambda lambda 0 0 0;...
    lambda lambda+2*mu lambda 0 0 0;...
    lambda lambda lambda+2*mu 0 0 0;...
    0 0 0 mu 0 0;...
    0 0 0 0 mu 0;...
    0 0 0 0 0 mu];
I = eye(3);
selectMapVoigt(mapVoigtObject,dimension,'symmetric');
% initialization energy and storageFEObject
strainEnergy = 0;
initializeDataFE(storageFEObject);
%% element loop for residual and tangent
for e = 1:numberOfElements
    Re = zeros(numberOfDofs,1);
    Ke = zeros(numberOfDofs);
    edN1 = qN1(edof(e,:),1:dimension)';
    J = qR(edof(e,:),1:dimension)'*dNr';
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
%         [Re1, Ke1] = gaussDisplacementSCSaintVenantEndpoint(edN1(:),k,dimension,J,dNrAll,DMat,gaussWeight,I);
        [Ke1,Re1] = admDiffFor(@gaussDisplacementSCSaintVenantEndpoint,1,edN1(:),k,dimension,J,dNr,DMat,gaussWeight,I,admOptions('flags', '--check-certificate'));
        Re = Re + Re1;
        Ke = Ke + Ke1(1:numberOfDofs,1:numberOfDofs); 
    end
    storeDataFE(storageFEObject,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end