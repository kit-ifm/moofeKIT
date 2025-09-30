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
solidObject.elementDisplacementType = 'eas';
solidObject.shapeFunctionObject.order = 1;
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshCooksMembrane(2, 2, solidObject.shapeFunctionObject.order, false);
solidObject.materialObject.name = 'HookeESZ';
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1;
solidObject.materialObject.nu = 1/3;
solidObject.dimension = 2;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order+1)^2;
solidObject.mixedFEObject.condensation = true;
solidObject.mixedFEObject.typeShapeFunctionData = 4;

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
dirichletObject.nodalDof = 1:2;

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [0; 1/16];
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)

%% postprocessing
% nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-48) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-60) < 1e-8); % upper right node
nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-48) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-52) < 1e-8); % middle right node
yDisplacementUpperRightNode = solidObject.qN1(nodeToMeasure,2)-solidObject.qR(nodeToMeasure,2);
fprintf('\nDisplacement: %4.3f\n',yDisplacementUpperRightNode)