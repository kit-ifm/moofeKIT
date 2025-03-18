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
lengthX = 1;
lengthY = 1;
axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
axisymmetricSolidObject.dimension = 2;

% axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
% axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
% axisymmetricSolidObject.numericalTangentObject.type = 'complex';

order = 2;
serendipity = true;
numberOfElementsX = 2;
numberOfElementsY = 2;
innerRadius = 1;
axisymmetricSolidObject.shapeFunctionObject.order = order;
axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

[axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
nodeToDistort = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;

if strcmp(meshType, 'B')
    axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) + 0.12;
    axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) + 0.11;
    axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);
end

% axisymmetricSolidObject.elementDisplacementType = 'eas';
% axisymmetricSolidObject.elementNameAdditionalSpecification = 'KasperTaylor';
% axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 2;
% axisymmetricSolidObject.mixedFEObject.condensation = true;

axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';

axisymmetricSolidObject.materialObject.name = 'Hooke';
axisymmetricSolidObject.materialObject.rho = 0;
axisymmetricSolidObject.materialObject.E = 250;
axisymmetricSolidObject.materialObject.nu = 0.1;

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
dirichletBoundary.nodalDof = 2;
dirichletBoundary.timeFunction = @(t, XYZ) -lengthY/2;
dirichletBoundary.masterObject = axisymmetricSolidObject;

dirichletBoundary2 = dirichletClass(dofObject);
dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
dirichletBoundary2.nodalDof = 2;
dirichletBoundary2.timeFunction = @(t, XYZ) lengthY/2;
dirichletBoundary2.masterObject = axisymmetricSolidObject;


% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = axisymmetricSolidObject;
neumannObject.loadVector = [1; 0]; % Komponenten des Kraftvektors anpassen
neumannObject.meshObject.edof = bounEdof.SX1;

%% solver
dofObject = runNewton(setupObject, dofObject);
nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 1) - axisymmetricSolidObject.qN(nodeToMeasure, 1);
disp(['End displacement: ', num2str(endDisplacement)]);
exactDisplacement = 1/750*(5+3*axisymmetricSolidObject.materialObject.nu-2*axisymmetricSolidObject.materialObject.nu^2);
disp(['Ratio FEM/Exact: ', num2str(endDisplacement/exactDisplacement)]);
