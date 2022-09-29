run('../prepareWorkspace');

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.fileName = 'patchTest';
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.plotFlag = true;
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
part1 = solidClass(dofObject);
[part1.nodes,part1.edof] = meshGeneratorCube(1,1,1,1,1,1,1,false);
part1.nodes = part1.nodes + 0.5;
% part1.materialName = 'Hooke'; 
% part1.materialName = 'SaintVenant';
part1.elementDisplacementType = 'displacementSC';
part1.materialName = 'MooneyRivlin';
part1.materialData.rho = 0;
part1.materialData.a = 1;
part1.materialData.b = 1;
part1.materialData.c = 1;
part1.materialData.d = 2*(part1.materialData.a + 2*part1.materialData.b);
part1.dimension = 3;
part1.orderShapeFunctions = 1;
part1.numberOfGausspoints = 8;

part3=solidClass(dofObject);

part2 = solidClass(dofObject);
[part2.nodes,part2.edof] = meshGeneratorCube(1,1,1,3,3,3,1,false);
part2.nodes = part2.nodes + 2.5;
part2.elementDisplacementType = 'displacementSC';
part2.materialName = 'MooneyRivlin';
part2.materialData.rho = 1;
part2.materialData.a = 1;
part2.materialData.b = 1;
part2.materialData.c = 1;
part2.materialData.d = 2*(part2.materialData.a + 2*part2.materialData.b);
part2.dimension = 3;
part2.orderShapeFunctions = 1;
part2.numberOfGausspoints = 8;

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

boundaryK2 = dirichletClass(dofObject);
boundaryK2.nodeList = find(part2.nodes(:,3) == 3);
boundaryK2.nodalDof = 3;
boundaryK2.masterObject = part2;
boundaryK2.timeFunction = str2func('@(t,Z) (Z - 0.25*t).*(t >= 0)');

boundaryK2(2) = dirichletClass(dofObject);
boundaryK2(2).nodeList = find(part2.nodes(:,3) == 2);
boundaryK2(2).nodalDof = 1:3;
boundaryK2(2).masterObject = part2;

%% solver
dofObject = runNewton(setupObject,dofObject);
