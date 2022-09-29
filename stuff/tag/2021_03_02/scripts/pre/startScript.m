run('../prepareWorkspace');

%% setup
setupObject(1) = setupClass;
setupObject(1).fileName = 'patchTest';
setupObject(1).totalTimeSteps = 4;

%% mainObjects
solidObject(1) = solidClass;
% TODO baseFEClass 
% TODO solidSuperClass (derive from baseFEClass); 

% meshGenerator
[solidObject(1).nodes,solidObject(1).edof] = meshGeneratorCube(1,1,1,20,20,20,1,false);
solidObject(1).nodes = solidObject(1).nodes + 0.5;
% TODO abaqusInterface

% material
solidObject(1).materialName = 'Hooke'; 
solidObject(1).materialData.rho = 0;
solidObject(1).materialData.lambda = 27.7777;
solidObject(1).materialData.mu = 41.6666;
solidObject(1).dimension = 3;
solidObject(1).orderShapeFunctions = 1;
solidObject(1).numberOfGausspoints = 8;

% % test
% solidObject(2) = solidClass;
% [solidObject(2).nodes,solidObject(2).edof] = meshGeneratorCube(1,1,1,3,3,3,1,false);
% solidObject(2).nodes = solidObject(2).nodes + 2.5;
% solidObject(2).materialName = 'Hooke';
% solidObject(2).materialData.rho = 1000;
% solidObject(2).materialData.lambda = 27.7777;
% solidObject(2).materialData.mu = 41.6666;
% solidObject(2).dimension = 3;
% solidObject(2).orderShapeFunctions = 1;
% solidObject(2).numberOfGausspoints = 8;

dirichletObject(1) = dirichletClass;
dirichletObject(1).masterNodesDofDirichlet = find(solidObject(1).nodes(:,3) == 0);
dirichletObject(1).masterNodesDirection = 3;
dirichletObject(1).masterObject = solidObject(1);

dirichletObject(2) = dirichletClass;
dirichletObject(2).masterNodesDofDirichlet = find(solidObject(1).nodes(:,1) == 0);
dirichletObject(2).masterNodesDirection = 1;
dirichletObject(2).masterObject = solidObject(1);

dirichletObject(3) = dirichletClass;
dirichletObject(3).masterNodesDofDirichlet = find(solidObject(1).nodes(:,2) == 0);
dirichletObject(3).masterNodesDirection = 2;
dirichletObject(3).masterObject = solidObject(1);

dirichletObject(4) = dirichletClass;
dirichletObject(4).masterNodesDofDirichlet = find(solidObject(1).nodes(:,3) == 1);
dirichletObject(4).masterNodesDirection = 3;
dirichletObject(4).masterObject = solidObject(1);
dirichletObject(4).timeFunction = @(t) (-0.5*t)*(t <= 1) + (-0.5)*(t > 1);

%% solver
tic
runNewton 
% TODO: evtl. als funktion oder gar als Object
toc