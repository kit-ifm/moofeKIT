%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
% setupObject.plotObject.time = 'R';
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 23;
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
orderShapeFunctions = 1;
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes, plateObject.meshObject.edof, edofBoundary] = meshCircularPlate(5, 3*4*4, orderShapeFunctions, false);
plateObject.materialObject.name = 'Hooke';
plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin';
% plateObject.elementNameAdditionalSpecification = 'SimoRifai';
% plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
% plateObject.elementDisplacementType = 'eas';
plateObject.mixedFEObject.typeShapeFunctionData = 4;
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 10.92e5;
plateObject.materialObject.nu = 0.3;
plateObject.h = 0.1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions + 1)^2;
% plateObject.elementDisplacementType = 'selectiveReducedIntegration';
% plateObject.mixedFEObject.condensation = true;

% Dirichlet boundary
numberOfNodes = size(plateObject.meshObject.nodes, 1);
nodeListDirichletBoundary = [];
for i = 1:numberOfNodes
    if abs(sqrt(plateObject.meshObject.nodes(i, 1)^2+plateObject.meshObject.nodes(i, 2)^2)-5) < 1e-10
        nodeListDirichletBoundary = [nodeListDirichletBoundary, i];
    end
end

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = nodeListDirichletBoundary;
boundary1.nodalDof = 1:3;
boundary1.masterObject = plateObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
boundary2.nodalDof = 2;
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
boundary3.nodeList = find(plateObject.meshObject.nodes(:, 2) == 0);
boundary3.nodalDof = 3;
boundary3.masterObject = plateObject;

% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = plateObject;
neumannObject.loadVector = -[5 * 2.048; 0; 0]; %-[2.511698113207548e-05; 0; 0]; % Komponenten des Kraftvektors anpassen
neumannObject.meshObject.edof = plateObject.meshObject.edof;

%% solver
dofObject = runNewton(setupObject, dofObject);
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(1, 1))))]);
