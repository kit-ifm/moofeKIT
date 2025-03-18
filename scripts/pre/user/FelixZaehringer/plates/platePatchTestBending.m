%% Plate Patch Test Bending
% expected result: 0
%
% REFERENCE
% https://doi.org/10.1108/eb023593

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'platePatchTestBending';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type', 'Cauchy', 'component', 23);
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
plateLength = 10;
orderShapeFunctions = 1;
serendipity = false;
plateObject = plateClass(dofObject);

[plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshPatchTestDistorted2D(plateLength, plateLength, orderShapeFunctions, serendipity);

plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 2.1e6;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1; % 0.001
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions + 1)^2;
plateObject.elementDisplacementType = 'eas'; % displacement, eas, selectiveReducedIntegration
plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin'; % PetrovGalerkinBatheDvorkin, BatheDvorkin, SimoRifai
% plateObject.elementDisplacementType = 'displacement'; % displacement, eas, selectiveReducedIntegration
% plateObject.elementNameAdditionalSpecification = 'BatheDvorkin'; % PetrovGalerkinBatheDvorkin, BatheDvorkin, SimoRifai
plateObject.mixedFEObject.condensation = false;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
boundary1.nodalDof = [1, 2];
boundary1.masterObject = plateObject;

% Load
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = plateObject;
nodalLoadObject.loadVector = [0; 10; 0]/2;
nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength);

%% solver
dofObject = runNewton(setupObject, dofObject);
