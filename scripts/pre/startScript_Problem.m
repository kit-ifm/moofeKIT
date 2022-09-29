%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 2;
setupObject.plotFlag = false;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
[solidObject.nodes,solidObject.edof] = meshGeneratorCube(1,1,1,3,3,3,2,true);
solidObject.nodes = solidObject.nodes + 0.5;
% solidObject.materialObject.name = 'SCNeoHooke';
% solidObject.materialObject.name = 'SaintVenant'; 
solidObject.materialObject.name = 'MooneyRivlin';
solidObject.elementDisplacementType = 'mixedSC';
% solidObject.elementDisplacementType = 'displacementSC';
% solidObject.elementDisplacementType = 'eas';
solidObject.materialObject.rho = 0;
% solidObject.materialObject.lambda = 27.7777;
% solidObject.materialObject.mu = 41.6666;
solidObject.materialObject.a = 1;
solidObject.materialObject.b = 1;
solidObject.materialObject.c = 1;
solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
solidObject.dimension = 3;
solidObject.orderShapeFunctions = 2;
solidObject.numberOfGausspoints = 8;
solidObject.mixedFEObject.condensation = true;
solidObject.mixedFEObject.orderShapeFunction = 1;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidObject.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = solidObject;

boundary1(2) = dirichletClass(dofObject);
boundary1(2).nodeList = find(solidObject.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = solidObject;

boundary1(3) = dirichletClass(dofObject);
boundary1(3).nodeList = find(solidObject.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = solidObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(solidObject.nodes(:,3) == 1);
boundary2.nodalDof = 3;
boundary2.masterObject = solidObject;
boundary2.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');
%% solver
tic
dofObject = runNewton(setupObject,dofObject);
plot(solidObject)
toc

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
strainEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
figure; 
plot(timeVector,kineticEnergy + strainEnergy);