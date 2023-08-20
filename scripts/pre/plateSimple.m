%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLength = 100;
orderShapeFunctions = 1;
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLength, plateLength, 4, 4, orderShapeFunctions);
plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength/2;
plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e5;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions+1)^2;
% plateObject.elementDisplacementType = 'selectiveReducedIntegration';
% plateObject.mixedFEObject.condensation = true;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

% boundary2 = dirichletClass(dofObject);
% boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
% boundary2.nodalDof = 2;
% boundary2.masterObject = plateObject;
% 
% boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
% boundary3.nodalDof = 3;
% boundary3.masterObject = plateObject;

% boundary1 = dirichletClass(dofObject);
% boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
% boundary1.nodalDof = 1;
% boundary1.masterObject = plateObject;
% 
% boundary2 = dirichletClass(dofObject);
% boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
% boundary2.nodalDof = 2;
% boundary2.masterObject = plateObject;
% 
% boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
% boundary3.nodalDof = 3;
% boundary3.masterObject = plateObject;

% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = plateObject;
neumannObject.loadVector = -[1; 0; 0];                         % Komponenten des Kraftvektors anpassen
neumannObject.loadGeometry = 'area';
neumannObject.meshObject.edof = plateObject.meshObject.edof;

% % Nodal Forces
% nodalForceObject = nodalForceClass(dofObject);
% nodalForceObject.masterObject = plateObject;
% nodalForceObject.loadVector = -[0; 0.5; 0];                         % Komponenten des Kraftvektors anpassen
% nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength));


%% solver
dofObject = runNewton(setupObject,dofObject);
%plot(plateObject, setupObject)
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);