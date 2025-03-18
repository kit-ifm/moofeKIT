%% MOMENT TEST
% REFERENCE: CTWM Script (W. Wagner)

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 11;
setupObject.usePreconditioning = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLengthX = 5;
plateLengthY = 5;
plateObject = plateClass(dofObject);
order = 1;
[plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, 5, 5, order, false);
plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;
plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 10.92e5;
plateObject.materialObject.nu = 0.3;
plateObject.h = 0.1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = order;
plateObject.shapeFunctionObject.numberOfGausspoints = 4;
plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
% plateObject.elementDisplacementType = 'eas';
% plateObject.elementDisplacementType = 'selectiveReducedIntegration';
% plateObject.mixedFEObject.condensation = true;
% plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0);
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX);
boundary2.nodalDof = 2;
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
boundary3.nodeList = find(plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary3.nodalDof = 3;
boundary3.masterObject = plateObject;

% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = plateObject;
neumannObject.loadGeometry = 'area';
neumannObject.loadVector = [-1; 0; 0];
neumannObject.meshObject.edof = plateObject.meshObject.edof;


%% solver
dofObject = runNewton(setupObject,dofObject);
%disp(['Central displacement: ', num2str(plateObject.qN1(find(plateObject.meshObject.nodes(:,1) == plateLengthX/2 & plateObject.meshObject.nodes(:,2) == plateLengthY/2), 1))]);
disp(['Central displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);
% Ziel EAS: 4.327875 (4x4), 4.1349466 (8x8), 40971.8884
%plot(plateObject, setupObject)