%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
% setupObject.plotObject.time = 'N';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLength = 100;
order = 1;
serendipity = false;

plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof, nodeListDirichlet] = meshRhombicPlate(plateLength, 4, order, serendipity);
plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 10^5;
plateObject.materialObject.nu = 0.3;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = order;
plateObject.shapeFunctionObject.numberOfGausspoints = (order+1)^2;
plateObject.elementDisplacementType = 'displacement'; % eas, displacement, easPetrovGalerkin, displacementPetrovGalerkin
plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
% plateObject.elementDisplacementType = 'selectiveReducedIntegration';
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;
% plateObject.numericalTangentObject.computeNumericalTangent = true;

% Dirichlet boundary
% boundary1 = dirichletClass(dofObject);
% boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == 100 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == 100);
% boundary1.nodalDof = 1;
% boundary1.masterObject = plateObject;
% 
% boundary2 = dirichletClass(dofObject);
% boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == 100 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == 100);
% boundary2.nodalDof = 2;
% boundary2.masterObject = plateObject;
% 
% boundary3 = dirichletClass(dofObject);
% boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,1) == 100 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == 100);
% boundary3.nodalDof = 3;
% boundary3.masterObject = plateObject;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = nodeListDirichlet;
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;

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
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = plateObject;
neumannObject.loadVector = [-1; 0; 0];                         % Komponenten des Kraftvektors anpassen
neumannObject.meshObject.edof = plateObject.meshObject.edof;

% % Nodal Forces
% nodalForceObject = nodalForceClass(dofObject);
% nodalForceObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
% nodalForceObject.masterObject = plateObject;
% nodalForceObject.forceVector = -[0; 0.5; 0];                         % Komponenten des Kraftvektors anpassen
% nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength));


%% solver
dofObject = runNewton(setupObject,dofObject);
disp(['Central displacement: ', num2str(plateObject.qN1(find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == 0), 1))]);
%plot(plateObject, setupObject)
% Ziel EAS: 3.93046 (4x4), 3.98412 (8x8)


function [nodes, edof, nodeListDirichletBoundary] = meshRhombicPlate(sideLength, numberOfElementsPerSide, order, serendipity)
%MESHRHOMBICPLATE Mesh for a rhombic plate.
%
%   REFERENCES
%   https://doi.org/10.1002/nme.1620290802, p.1629
%
%   CREATOR(S)
%   Felix Zaehringer

angle = 30; % angle of inclination in degrees

lengthX = sideLength;
lengthY = sind(angle) * sideLength;
[nodes, edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsPerSide, numberOfElementsPerSide, order, serendipity);

nodeListDirichletBoundary = find(nodes(:, 1) == -lengthX/2 | nodes(:, 1) == lengthX/2 | nodes(:, 2) == -lengthY/2 | nodes(:, 2) == lengthY/2);

nodes(:, 1) = nodes(:, 1) + tand(90-angle) * nodes(:, 2);

end
