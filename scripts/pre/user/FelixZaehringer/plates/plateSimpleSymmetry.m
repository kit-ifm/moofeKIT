%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.usePreconditioning = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLengthX = 50;
plateLengthY = 50;
plateObject = plateClass(dofObject);
order = 1;
[plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, 50, 50, order);
plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;
plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e4;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = order;
plateObject.shapeFunctionObject.numberOfGausspoints = 4;
% plateObject.elementDisplacementType = 'eas';
% plateObject.elementDisplacementType = 'selectiveReducedIntegration';
plateObject.mixedFEObject.condensation = true;
% plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0);
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

boundary2 = dirichletClass(dofObject);
% boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,1) == plateLengthX);
boundary2.nodeList = find(plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,1) == plateLengthX);
boundary2.nodalDof = 2;
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary3.nodalDof = 3;
boundary3.masterObject = plateObject;

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

% % % Neumann boundary
% neumannObject = neumannClass(dofObject);
% neumannObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
% neumannObject.masterObject = plateObject;
% neumannObject.forceVector = [-1; 0; 0];                         % Komponenten des Kraftvektors anpassen
% neumannObject.shapeFunctionObject.dimension = plateObject.dimension;
% neumannObject.shapeFunctionObject.order = plateObject.shapeFunctionObject.order;
% neumannObject.shapeFunctionObject.numberOfGausspoints = plateObject.shapeFunctionObject.numberOfGausspoints;  % Anzahl der Gaußpunkten auf dem Neumannrand anpassen
% neumannObject.meshObject.edof = plateObject.meshObject.edof;

% Nodal Forces
nodalForceObject = nodalForceClass(dofObject);
nodalForceObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
nodalForceObject.masterObject = plateObject;
nodalForceObject.forceVector = -[16.367; 0; 0];                         % Komponenten des Kraftvektors anpassen
% nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthX));
nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX & plateObject.meshObject.nodes(:,2) == plateLengthY);

%% solver
dofObject = runNewton(setupObject,dofObject);
%disp(['Central displacement: ', num2str(plateObject.qN1(find(plateObject.meshObject.nodes(:,1) == plateLengthX/2 & plateObject.meshObject.nodes(:,2) == plateLengthY/2), 1))]);
disp(['Central displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);
% Ziel EAS: 4.327875 (4x4), 4.1349466 (8x8), 40971.8884
%plot(plateObject, setupObject)