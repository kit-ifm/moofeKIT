run('../prepareWorkspace');

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.fileName = 'patchTest';
setupObject.totalTimeSteps = 4;
setupObject.totalTime = 1;
setupObject.plotFlag = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
part1 = solidClass(dofObject);
[part1.nodes,part1.edof] = meshGeneratorCube(1,1,1,2,2,2,1,false);
part1.nodes = part1.nodes + 0.5;
% part1.materialName = 'Hooke'; 
part1.materialName = 'SaintVenant';
part1.materialData.rho = 0;
part1.materialData.lambda = 27.7777;
part1.materialData.mu = 41.6666;
part1.dimension = 3;
part1.orderShapeFunctions = 1;
part1.numberOfGausspoints = 8;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(part1.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = part1;

boundary1(2) = dirichletClass(dofObject);
boundary1(2).nodeList = find(part1.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = part1;

boundary1(3) = dirichletClass(dofObject);
boundary1(3).nodeList = find(part1.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = part1;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(part1.nodes(:,3) == 1);
boundary2.nodalDof = 3;
boundary2.masterObject = part1;
boundary2.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

%% solver
dofObject = runNewton(setupObject,dofObject);
