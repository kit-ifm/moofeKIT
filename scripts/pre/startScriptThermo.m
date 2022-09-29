%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'temp';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
thermoObject = thermoClass(dofObject);
length = 0.5;
orderShapeFunctions = 1;
[thermoObject.meshObject.nodes,thermoObject.meshObject.edof] = meshOneDimensional(length, 3, orderShapeFunctions);
thermoObject.meshObject.nodes = thermoObject.meshObject.nodes + length/2;
thermoObject.elementDisplacementType = 'thermo';
thermoObject.materialObject.name = 'LinearFourier';
setupObject.integrator = 'Endpoint';
thermoObject.materialObject.rho = 0;
thermoObject.dimension = 1;
thermoObject.shapeFunctionObject.order = orderShapeFunctions;
thermoObject.shapeFunctionObject.numberOfGausspoints = 1;

% dirichlet boundaries
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(thermoObject.meshObject.nodes(:,1) == 0);
dirichletObject1.nodalDof = 1;
dirichletObject1.timeFunction = @(t, Z) -20;
dirichletObject1.masterObject = thermoObject;

dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = find(thermoObject.meshObject.nodes(:,1) == 0.5);
dirichletObject2.nodalDof = 1;
dirichletObject2.timeFunction = @(t, Z) 30;
dirichletObject2.masterObject = thermoObject;

%% solver
dofObject = runNewton(setupObject,dofObject);

%% plot temperature
x = thermoObject.meshObject.nodes;
y = thermoObject.qN1;
figure;
plot(x,y);