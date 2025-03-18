%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'continuum2DTestDistorted';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 33;
setupObject.plotObject.view = 2;
setupObject.plotObject.time = 'R';
setupObject.integrator = 'Endpoint';
setupObject.newton.tolerance = 1e-4;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
meshType = 'A'; % A or B (regular/distorted)
lengthX = 6;
lengthY = 1;
axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
axisymmetricSolidObject.dimension = 2;

% axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
% axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
% axisymmetricSolidObject.numericalTangentObject.type = 'complex';

order = 1;
serendipity = false;
numberOfElementsX = 5;
numberOfElementsY = 1;
innerRadius = 3;
axisymmetricSolidObject.shapeFunctionObject.order = order;
axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

[axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;
axisymmetricSolidObject.meshObject.nodes([2, 2+6], 1) = innerRadius + 0.5;
axisymmetricSolidObject.meshObject.nodes([3, 3+6], 1) = innerRadius + 5/4;
axisymmetricSolidObject.meshObject.nodes([4, 4+6], 1) = innerRadius + 9/4;
axisymmetricSolidObject.meshObject.nodes([5, 5+6], 1) = innerRadius + 15/4;

if strcmp(meshType, 'B')
    axisymmetricSolidObject.meshObject.nodes(2, 1) = axisymmetricSolidObject.meshObject.nodes(2, 1) - 1/6;
    axisymmetricSolidObject.meshObject.nodes(2+6, 1) = axisymmetricSolidObject.meshObject.nodes(2+6, 1) + 1/6;
    axisymmetricSolidObject.meshObject.nodes(3, 1) = axisymmetricSolidObject.meshObject.nodes(3, 1) + 1/6;
    axisymmetricSolidObject.meshObject.nodes(3+6, 1) = axisymmetricSolidObject.meshObject.nodes(3+6, 1) - 1/6;
    axisymmetricSolidObject.meshObject.nodes(4, 1) = axisymmetricSolidObject.meshObject.nodes(4, 1) + 1/3;
    axisymmetricSolidObject.meshObject.nodes(4+6, 1) = axisymmetricSolidObject.meshObject.nodes(4+6, 1) - 1/3;
    axisymmetricSolidObject.meshObject.nodes(5, 1) = axisymmetricSolidObject.meshObject.nodes(5, 1) - 1/3;
    axisymmetricSolidObject.meshObject.nodes(5+6, 1) = axisymmetricSolidObject.meshObject.nodes(5+6, 1) + 1/3;
end
axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);

axisymmetricSolidObject.elementDisplacementType = 'eas';
axisymmetricSolidObject.elementNameAdditionalSpecification = 'SimoRifai';
axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 5;
axisymmetricSolidObject.mixedFEObject.condensation = true;

% axisymmetricSolidObject.elementDisplacementType = 'eas';
% axisymmetricSolidObject.elementNameAdditionalSpecification = 'KasperTaylor';
% axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 2;
% axisymmetricSolidObject.mixedFEObject.condensation = true;

% axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';

axisymmetricSolidObject.materialObject.name = 'Hooke';
axisymmetricSolidObject.materialObject.rho = 0;
axisymmetricSolidObject.materialObject.E = 1000;
axisymmetricSolidObject.materialObject.nu = 0.25;

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
exactDisplacement = (1+axisymmetricSolidObject.materialObject.nu)*1*innerRadius^2/(axisymmetricSolidObject.materialObject.E*((innerRadius+lengthX)^2-innerRadius^2))*((innerRadius+lengthX)^2/innerRadius+(1-2*axisymmetricSolidObject.materialObject.nu)*innerRadius);
disp(['Ratio FEM/Exact: ', num2str(endDisplacement/exactDisplacement)]);
