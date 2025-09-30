%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembrane';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 3;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
% solidObject.elementDisplacementType = 'displacement';
solidObject.elementDisplacementType = 'eas';
solidObject.shapeFunctionObject.order = 1;
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshCooksMembrane(2, 2, solidObject.shapeFunctionObject.order, false);
solidObject.materialObject.name = 'HookeESZ';
% solidObject.materialObject.name = 'Hooke';
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1e2;
solidObject.materialObject.nu = 0.499;
solidObject.dimension = 2;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;
solidObject.mixedFEObject.typeShapeFunctionData = 4;
solidObject.mixedFEObject.condensation = false;

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = 1:2;
dirichletObject.timeFunction = str2func('@(t) 0');

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [0; 5];
neumannObject.timeFunction = @(t) t;
neumannObject.meshObject.edof = edofBoundary;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)