%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'continuum2DTestDistorted';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 33;
setupObject.plotObject.view = 2;
setupObject.integrator = 'Endpoint';
setupObject.newton.tolerance = 1e-4;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
meshType = 'B'; % A or B (regular/distorted)
lengthX = 5;
lengthY = 0.1; % 0.1
axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
axisymmetricSolidObject.dimension = 2;

% axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
% axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
% axisymmetricSolidObject.numericalTangentObject.type = 'complex';

order = 2;
serendipity = true;
numberOfElementsX = 8;
numberOfElementsY = 1;
s = lengthX/numberOfElementsX*0.4;
innerRadius = 0;
axisymmetricSolidObject.shapeFunctionObject.order = order;
axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

if strcmp(meshType, 'A')
    [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
    axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;
elseif strcmp(meshType, 'B')
    [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
    axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;
    nodesToDistortBottom = find(axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    nodesToDistortTop = find(axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
    if order == 2
        nodesToDistortBottom = nodesToDistortBottom(2:2:end);
        nodesToDistortTop = nodesToDistortTop(2:2:end);
    end
    axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) + s;
    axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) - s;
%     axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) - s;
%     axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) + s;
elseif strcmp(meshType, 'C')
    [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity);
    axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + innerRadius;
    axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) - lengthY/2;
end
axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);

axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
axisymmetricSolidObject.materialObject.name = 'Hooke';
axisymmetricSolidObject.materialObject.rho = 0;
axisymmetricSolidObject.materialObject.E = 10.92e5;
axisymmetricSolidObject.materialObject.nu = 0.3;

% Dirichlet boundary
boundaryCondition = 'clamped'; % clamped or simplySupported or doubleSupport
if strcmp(boundaryCondition, 'clamped')
    dirichletBoundary = dirichletClass(dofObject);
    dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX);
    dirichletBoundary.nodalDof = 1;
    dirichletBoundary.timeFunction = @(t, XYZ) lengthX;
    dirichletBoundary.masterObject = axisymmetricSolidObject;
    
    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    dirichletBoundary2.nodalDof = 2;
    dirichletBoundary2.timeFunction = @(t, XYZ) -lengthY/2;
    dirichletBoundary2.masterObject = axisymmetricSolidObject;
    
    dirichletBoundary3 = dirichletClass(dofObject);
    dirichletBoundary3.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
    dirichletBoundary3.nodalDof = 2;
    dirichletBoundary3.timeFunction = @(t, XYZ) lengthY/2;
    dirichletBoundary3.masterObject = axisymmetricSolidObject;
elseif strcmp(boundaryCondition, 'simplySupported')
    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    dirichletBoundary2.nodalDof = 2;
    dirichletBoundary2.timeFunction = @(t, XYZ) -lengthY/2;
    dirichletBoundary2.masterObject = axisymmetricSolidObject;
elseif strcmp(boundaryCondition, 'doubleSupport')
    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    dirichletBoundary2.nodalDof = 2;
    dirichletBoundary2.timeFunction = @(t, XYZ) -lengthY/2;
    dirichletBoundary2.masterObject = axisymmetricSolidObject;

    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    dirichletBoundary2.nodalDof = 1;
    dirichletBoundary2.timeFunction = @(t, XYZ) lengthX;
    dirichletBoundary2.masterObject = axisymmetricSolidObject;
end

% Neumann boundary
loadCondition = 'lineLoad'; % lineLoad or endMoment
if strcmp(loadCondition, 'lineLoad')
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = axisymmetricSolidObject;
    neumannObject.loadVector = [0; -10.24]; % Komponenten des Kraftvektors anpassen
    neumannObject.meshObject.edof = bounEdof.SY2;
elseif strcmp(loadCondition, 'endMoment')
    M0 = 100;
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = axisymmetricSolidObject;
    nodalLoadObject.loadVector = M0/lengthY*[-1; 0]*2*pi;
    nodalLoadObject.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:,1) == innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:,2) == lengthY/2);
    
    nodalLoadObject2 = nodalLoadClass(dofObject);
    nodalLoadObject2.masterObject = axisymmetricSolidObject;
    nodalLoadObject2.loadVector = -M0/lengthX*[1; 0]*2*pi;
    nodalLoadObject2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:,1) == innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:,2) == -lengthY/2);
end

%% solver
dofObject = runNewton(setupObject, dofObject);
nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 2) - axisymmetricSolidObject.qN(nodeToMeasure, 2);
disp(['End displacement: ', num2str(endDisplacement)]);
