run('../prepareWorkspace');

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.fileName = 'solidThermo';
setupObject.totalTimeSteps = 4;
setupObject.totalTime = 1;
setupObject.plotFlag = true;

%% continuum Objects
part2 = solidThermoClass;
[part2.nodes,part2.edof] = meshGeneratorCube(1,1,1,2,2,2,1,false);
part2.nodes = part2.nodes + 2.5;
part2.nodes = [part2.nodes, 293.15*ones(size(part2.nodes,1),1)];
part2.elementDisplacementType = 'displacementSC';
part2.materialName = 'MooneyRivlin';
part2.materialData.rho = 1;
part2.materialData.a = 1;
part2.materialData.b = 1;
part2.materialData.c = 1;
part2.materialData.d = 2*(part2.materialData.a + 2*part2.materialData.b);
part2.materialData.kappa = 1;
part2.materialData.beta = 2.233*10^(-4);    % coupling parameter
part2.materialData.theta0 = 293.15;         % reference Temperature
part2.materialData.k0 = 0.1;                % thermal conductivity
part2.dimension = 3;
part2.orderShapeFunctions = 1;
part2.numberOfGausspoints = 8;

boundaryK2 = dirichletClass;
boundaryK2.nodeList = find(part2.nodes(:,3) == 3);
boundaryK2.nodalDof = 3;
boundaryK2.masterObject = part2;
boundaryK2.timeFunction = str2func('@(t,Z) (Z - 0.25*t).*(t >= 0)');

%% solver
dofObject = dofClass;   % required object for dof and object handling
clearvars('-except','setupObject','dofObject',dofObject.listContinuumObjects{1:end})

dofObject = runNewton(setupObject,dofObject);
