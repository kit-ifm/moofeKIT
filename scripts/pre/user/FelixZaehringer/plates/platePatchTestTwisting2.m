%% Plate Patch Test Twisting II
% expected result: 0.5
%
% REFERENCE
% https://doi.org/10.1108/eb023593

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'platePatchTestTwistingII';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type', 'Cauchy', 'component',12);
setupObject.newton.tolerance = 1e-5;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
plateLength = 10;
orderShapeFunctions = 1;
plateObject = plateClass(dofObject);

% [plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshPatchTestDistorted2D(plateLength, plateLength, orderShapeFunctions, serendipity);
[plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshRectangle(plateLength, plateLength, 10, 10, orderShapeFunctions, serendipity);
plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength/2;

plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 2.1e6;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1e-2; % 0.001
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions + 1)^2;
plateObject.elementDisplacementType = 'eas'; % displacement, eas, selectiveReducedIntegration
plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin'; % PetrovGalerkinBatheDvorkin, BatheDvorkin, SimoRifai
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find((plateObject.meshObject.nodes(:, 1) == 0 & (plateObject.meshObject.nodes(:, 2) == 0 | plateObject.meshObject.nodes(:, 2) == plateLength)) | (plateObject.meshObject.nodes(:, 1) == plateLength & plateObject.meshObject.nodes(:, 2) == 0));
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

% Load
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = plateObject;
nodalLoadObject.loadVector = -[1; 0; 0];
nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength & plateObject.meshObject.nodes(:,2) == plateLength);

%% solver
dofObject = runNewton(setupObject, dofObject);
