%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
[solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,1,1,2,2,2,1,false);
solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
% solidObject.materialObject.name = 'NeoHooke';
solidObject.materialObject.name = 'SaintVenant'; 
% solidObject.materialObject.name = 'SaintVenantAutoDiff'; 
% solidObject.elementDisplacementType = 'displacement';
solidObject.elementDisplacementType = 'displacementSC';
% solidObject.elementDisplacementType = 'easSC';
solidObject.materialObject.rho = 1;
solidObject.materialObject.lambda = 27.7777;
solidObject.materialObject.mu = 41.6666;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;
solidObject.mixedFEObject.condensation = false;
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = false;
% solidObject.numericalTangentObject.type = 'complex';

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidObject.meshObject.nodes(:,3) == 1);
boundary1.nodalDof = 3;
boundary1.masterObject = solidObject;
boundary1.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

%% solver
tic
dofObject = runNewton(setupObject,dofObject);
plot(solidObject, setupObject)
toc
