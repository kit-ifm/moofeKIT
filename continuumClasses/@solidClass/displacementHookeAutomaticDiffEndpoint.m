function displacementHookeAutomaticDiffEndpoint(obj,setupObject)
%% setup
% load objects
meshObject = obj.meshObject;
lagrangeShapeFunctionObject = obj.lagrangeShapeFunctionObject;
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
gaussWeight = lagrangeShapeFunctionObject.gaussWeight;
numberOfGausspoints = lagrangeShapeFunctionObject.numberOfGausspoints;
%   N = lagrangeShapeFunctionObject.N;
dNr = lagrangeShapeFunctionObject.dNr;
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
for e = 1:numberOfElements
    Re = zeros(numberOfDofs,1);
    Ke = zeros(numberOfDofs);
    edN1 = qN1(edof(e,:),1:dimension)';
    edRef = qR(edof(e,:),1:dimension)';
    uN1 = edN1(:) - edRef(:);
    J = qR(edof(e,:),1:dimension)'*dNr';
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
%         [Re1,Ke1] = gaussDisplacementHookeEnpoint(uN1,k,dimension,J,dNrAll,DMat,gaussWeight);        
%         [Ke1 Re1] = admDiffRev(@gaussDisplacementHookeEnpoint,1,uN1,k,dimension,J,dNrAll,DMat,gaussWeight,admOptions('flags', '--check-certificate'));
%         [a_uN1 a_k a_dimension a_J a_dNrAll a_DMat a_gaussWeight nr_Re] = a_gaussDisplacementHookeEnpoint(uN1, k, dimension, J, dNrAll, DMat, gaussWeight, eye(numberOfDOFs));
        [Ke1,Re1] = admDiffFor(@gaussDisplacementHookeEnpoint,1,uN1,k,dimension,J,dNr,DMat,gaussWeight,admOptions('flags', '--check-certificate'));
%         [Ke1, Re1]= g_gaussDisplacementHookeEnpointMF(1, uN1, 0*k, k, 0*dimension, dimension, 0*J, J, 0*dNrAll, dNrAll, 0*DMat, DMat, 0*gaussWeight, gaussWeight);
        Re = Re + Re1;
        Ke = Ke + Ke1(:,1:numberOfDofs);
%         Ke = Ke + Ke1;%(:,1:numberOfDOFs);
    end 
    storeDataFE(storageFEObject,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end