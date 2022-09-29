%% changelog:
% core changes
% [x] auto updates
% [x] flags in classes and runNewton-script for update-procedures (idea: for further classes hopefully no need for adjustement of runNewton-script)
% [failed] catch objects via property of classes by creation of an instance
% [x] dirichletClass -> implement full features (arbritary dof, timeFunction: include nodal information)
% addional changes
% [] baseFEClass (edof; globa.lTotalEdof)
% [] solidSuperClass (derive from baseFEClass);
% [] speed up: save dofsObject -> object search
% [] destruct all 'pre' Objects except continuumObjects
% [] dofs -> dof (e.g. dofsObject -> dofObject)
% [] generate dofsObject via setupClass in the background
% [] push + tag

%% TODO:
% [] mixed finite elements -> solver/static condensation
% [] solidThermoClass (proof of concept for auto updates)

%% TODO later:
% [] universal numerical tangent (FD & complex tangent)
% [] define testroutine implement git-CI (matlab-runner)
% [] check function- (or class-)structure of runNewton-script
% [] implement GFE routines
% [] implement FEidF routines
% [] implement KEuG routines
% [] implement abaqus-input-file interface
% [] implement class for startScript config

run('../prepareWorkspace');

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.fileName = 'patchTest';
setupObject.totalTimeSteps = 4;
setupObject.plotFlag = true;

dofsObject = dofsClass;

%% continuum Objects
part1 = solidClass(dofsObject);
[part1.nodes,part1.edof] = meshGeneratorCube(1,1,1,2,2,2,1,false);
part1.nodes = part1.nodes + 0.5;

% material
part1.materialName = 'Hooke'; 
part1.materialData.rho = 0;
part1.materialData.lambda = 27.7777;
part1.materialData.mu = 41.6666;
part1.dimension = 3;
part1.orderShapeFunctions = 1;
part1.numberOfGausspoints = 8;

part2 = solidClass(dofsObject);
[part2.nodes,part2.edof] = meshGeneratorCube(1,1,1,3,3,3,1,false);
part2.nodes = part2.nodes + 2.5;
part2.materialName = 'SaintVenant';
part2.materialData.rho = 0;
part2.materialData.lambda = 27.7777;
part2.materialData.mu = 41.6666;
part2.dimension = 3;
part2.orderShapeFunctions = 1;
part2.numberOfGausspoints = 8;
 
boundary1 = dirichletClass;
boundary1.nodeList = find(part1.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = part1;

boundary1(2) = dirichletClass;
boundary1(2).nodeList = find(part1.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = part1;

boundary1(3) = dirichletClass;
boundary1(3).nodeList = find(part1.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = part1;

boundary2 = dirichletClass;
boundary2.nodeList = find(part1.nodes(:,3) == 1);
boundary2.nodalDof = 3;
boundary2.masterObject = part1;
boundary2.timeFunction = str2func('@(t,Z) (Z - 0.5*t).*(t >= 0)');

boundaryK2 = dirichletClass;
boundaryK2.nodeList = find(part2.nodes(:,3) == 3);
boundaryK2.nodalDof = 3;
boundaryK2.masterObject = part2;
boundaryK2.timeFunction = str2func('@(t,Z) (Z - 0.25*t).*(t >= 0)');

boundaryK2(2) = dirichletClass;
boundaryK2(2).nodeList = find(part2.nodes(:,3) == 2);
boundaryK2(2).nodalDof = 1:3;
boundaryK2(2).masterObject = part2;

%% solver
runNewton 