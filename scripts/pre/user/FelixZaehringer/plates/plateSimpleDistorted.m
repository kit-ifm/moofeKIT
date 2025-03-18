%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type','Cauchy', 'component', 11);
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLength = 50;
s=0;
orderShapeFunctions = 1;
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof, boundaryEdof] = meshRectangle(plateLength, plateLength, 2, 2, orderShapeFunctions);
nodeToDistort = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == 0);
plateObject.meshObject.nodes(nodeToDistort, :) = plateObject.meshObject.nodes(nodeToDistort, :) + s;
plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength/2;
plateObject.materialObject.name = 'HookeBatheDvorkin'; %HookeBatheDvorkin
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e5;
plateObject.materialObject.nu = 0;%0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions+1)^2;
plateObject.elementDisplacementType = 'displacementPetrovGalerkin'; % displacementPetrovGalerkin, eas, selectiveReducedIntegration
% plateObject.elementNameAdditionalSpecification = 'SimoRifai';
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
% boundary1 = dirichletClass(dofObject);
% boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
% boundary1.nodalDof = 1;
% boundary1.masterObject = plateObject;

% boundary2 = dirichletClass(dofObject);
% boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
% boundary2.nodalDof = 2;
% boundary2.masterObject = plateObject;
% 
% boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == plateLength | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength);
% boundary3.nodalDof = 3;
% boundary3.masterObject = plateObject;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
boundary1.nodalDof = 1:3;
boundary1.masterObject = plateObject;
% 
boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength);
boundary2.nodalDof = 2;
boundary2.timeFunction = @(t,z) 1;
boundary2.masterObject = plateObject;
% 
% boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
% boundary3.nodalDof = 3;
% boundary3.masterObject = plateObject;

% % Neumann boundary
% neumannObject = neumannClass(dofObject);
% neumannObject.loadGeometry = 'area';                                % Dimension des zu betrachtenden Körpers angeben
% neumannObject.masterObject = plateObject;
% neumannObject.loadVector = -[1; 0; 0];                         % Komponenten des Kraftvektors anpassen
% neumannObject.meshObject.edof = plateObject.meshObject.edof;

% neumannObject = neumannClass(dofObject);
% neumannObject.loadGeometry = 'line';                                % Dimension des zu betrachtenden Körpers angeben
% neumannObject.masterObject = plateObject;
% neumannObject.loadVector = -[0; 20; 0];                         % Komponenten des Kraftvektors anpassen
% neumannObject.meshObject.edof = boundaryEdof.SX2;

% neumannObject2 = neumannClass(dofObject);
% neumannObject2.loadGeometry = 'line';                                % Dimension des zu betrachtenden Körpers angeben
% neumannObject2.masterObject = plateObject;
% neumannObject2.loadVector = [0; 0; 10];                         % Komponenten des Kraftvektors anpassen
% neumannObject2.meshObject.edof = boundaryEdof.SY1;

% % Nodal Forces
% nodalForceObject = nodalForceClass(dofObject);
% nodalForceObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
% nodalForceObject.masterObject = plateObject;
% nodalForceObject.forceVector = -[0; 0.5; 0];                         % Komponenten des Kraftvektors anpassen
% nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength));


%% solver
dofObject = runNewton(setupObject,dofObject);
%plot(plateObject, setupObject)
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);