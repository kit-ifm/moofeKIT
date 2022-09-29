%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembrane';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.dimension = 2;
solidObject.elementDisplacementType = 'selectiveReducedIntegration';
solidObject.shapeFunctionObject.order = 1;
numberOfElements = 10;
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshCooksMembrane(numberOfElements, numberOfElements, solidObject.shapeFunctionObject.order);
solidObject.materialObject.name = 'HookeEVZ';
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1e2;
solidObject.materialObject.nu = 0.499;
solidObject.mixedFEObject.condensation = true;
solidObject.mixedFEObject.typeShapeFunctionData = 0;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order+1)^2;

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.dimension = 2;
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
dirichletObject.nodalDof = 1:2;
dirichletObject.timeFunction = str2func('@(t) 0');

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.dimension = 2;
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidObject;
neumannObject.forceVector = [0;5];
neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = neumannObject.shapeFunctionObject.order+1;
neumannObject.projectionType = 'none';
neumannObject.timeFunction = @(t) t;
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end,2)-solidObject.qR(end,2);
fprintf('\nDisplacement: %4.3f\n',yDisplacementUpperRightNode)