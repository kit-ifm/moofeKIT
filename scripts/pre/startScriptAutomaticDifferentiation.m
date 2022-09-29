addpath /home/marlon/Projekte/moofeKIT/adimat-0.6.6-5529-DebianBuster-x86_64/
ADiMat_startup

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTestAutomaticDifferentiation';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.plotFlag = false;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
% [solidObject.nodes,solidObject.edof] = meshGeneratorCube(1,1,1,6,6,6,1,false);
[solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,1,1,2,2,2,1,false);
solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
% solidObject.materialObject.name = 'HookeAutomaticDiff'; 
% solidObject.materialObject.name = 'Hooke'; 
solidObject.materialObject.name = 'SaintVenantAutomaticDiff';
% solidObject.materialObject.name = 'SaintVenant';
% solidObject.elementDisplacementType = 'displacement';
solidObject.elementDisplacementType = 'displacementSC';
solidObject.materialObject.rho = 0;
solidObject.materialObject.lambda = 27.7777;
solidObject.materialObject.mu = 41.6666;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.orderShapeFunctions = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;
%% FIXME automatisch anstupsen???
solidObject.shapeFunctionObject.computeLagrangeShapeFunction(size(solidObject.meshObject.edof,2),solidObject.dimension);

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidObject.meshObject.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = solidObject;

boundary1(2) = dirichletClass(dofObject);
boundary1(2).nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = solidObject;

boundary1(3) = dirichletClass(dofObject);
boundary1(3).nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = solidObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(solidObject.meshObject.nodes(:,3) == 1);
boundary2.nodalDof = 3;
boundary2.masterObject = solidObject;
boundary2.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

%% solver
dofObject = runNewton(setupObject,dofObject);

%% post
% filename = 'execute';
% VTKPlot(filename,'unstructured_grid',part1.qN1,part1.edof)
