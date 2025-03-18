%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-6;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLength = 100;
orderShapeFunctions = 1;
serendipity = true;
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLength, plateLength, 2, 2, orderShapeFunctions, serendipity);
plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength/2;
plateObject.meshObject.nodes(:, 1) = plateObject.meshObject.nodes(:, 1) + tand(90-30) * plateObject.meshObject.nodes(:, 2);
plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e5;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions+1)^2;
plateObject.elementDisplacementType = 'displacement'; % easPetrovGalerkin, eas, selectiveReducedIntegration, displacement
plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = plateObject;
neumannObject.loadGeometry = 'area';
neumannObject.loadVector = -[1; 0; 0];                         % Komponenten des Kraftvektors anpassen
neumannObject.meshObject.edof = plateObject.meshObject.edof;

%% solver
dofObject = runNewton(setupObject,dofObject);
%plot(plateObject, setupObject)
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);