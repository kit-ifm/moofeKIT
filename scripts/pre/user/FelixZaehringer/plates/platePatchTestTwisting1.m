%% Plate Patch Test Twisting I
% expected result: 1
%
% REFERENCE
% https://doi.org/10.1108/eb023593

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'platePatchTestTwistingI';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type', 'Cauchy', 'component', 12);
setupObject.newton.tolerance = 1e-5;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
plateLength = 10;
orderShapeFunctions = 1;
serendipity = false;
plateObject = plateClass(dofObject);

% [plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshPatchTestDistorted2D(plateLength, plateLength, orderShapeFunctions, serendipity);
[plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshRectangle(plateLength, plateLength, 1, 1, orderShapeFunctions, serendipity);
plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength/2;

plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 2.1e6;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions + 1)^2;
plateObject.elementDisplacementType = 'eas'; % displacement, eas, selectiveReducedIntegration
plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin'; % PetrovGalerkinBatheDvorkin, BatheDvorkin, SimoRifai
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0 | (plateObject.meshObject.nodes(:, 1) == plateLength & plateObject.meshObject.nodes(:, 2) == 0));
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
boundary2.nodalDof = 2;
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
boundary3.nodeList = find(plateObject.meshObject.nodes(:, 2) == 0);
boundary3.nodalDof = 3;
boundary3.masterObject = plateObject;

% Load
nodalLoadObject1 = nodalLoadClass(dofObject);
nodalLoadObject1.masterObject = plateObject;
nodalLoadObject1.loadVector = [0; 0; 10]/2*1000;
nodalLoadObject1.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength);

nodalLoadObject2 = nodalLoadClass(dofObject);
nodalLoadObject2.masterObject = plateObject;
nodalLoadObject2.loadVector = [0; 10; 0]/2*1000;
nodalLoadObject2.nodeList = find(plateObject.meshObject.nodes(:,2) == plateLength);

%% solver
dofObject = runNewton(setupObject, dofObject);
