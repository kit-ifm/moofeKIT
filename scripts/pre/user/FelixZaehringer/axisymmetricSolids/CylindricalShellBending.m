%% settings for simulation
M0 = 2000;
a = 167.5;
nu = 0.499999;
t = 1;
EModul = 11250;
L = 51;

beta4_a = 3*(1-nu^2)/(a^2*t^2);
beta4 = EModul*t/(4*a^2*D);
beta2 = sqrt(beta4);

D = EModul*t^3/(12*(1-nu^2));

uAnalytical = M0/(2*beta2*D);

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'continuum2DTestDistorted';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'zero';
setupObject.plotObject.stress.component = 11;
setupObject.plotObject.view = 2;
setupObject.integrator = 'Endpoint';
setupObject.newton.tolerance = 1e-4;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
meshType = 'A'; % A or B (regular/distorted)
lengthX = t;
lengthY = L;
axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
axisymmetricSolidObject.dimension = 2;

% axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
% axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
% axisymmetricSolidObject.numericalTangentObject.type = 'complex';

order = 1;
serendipity = false;
numberOfElementsX = 1;
numberOfElementsY = 17;
innerRadius = a-lengthX/2;
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

% axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';

axisymmetricSolidObject.materialObject.name = 'Hooke';
axisymmetricSolidObject.materialObject.rho = 0;
axisymmetricSolidObject.materialObject.E = EModul;
axisymmetricSolidObject.materialObject.nu = 0.1;

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
dirichletBoundary.nodalDof = 2;
dirichletBoundary.timeFunction = @(t, XYZ) -lengthY/2;
dirichletBoundary.masterObject = axisymmetricSolidObject;

dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2 & axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius);
dirichletBoundary.nodalDof = 1;
dirichletBoundary.timeFunction = @(t, XYZ) innerRadius;
dirichletBoundary.masterObject = axisymmetricSolidObject;

dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2 & axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX);
dirichletBoundary.nodalDof = 1;
dirichletBoundary.timeFunction = @(t, XYZ) innerRadius+lengthX;
dirichletBoundary.masterObject = axisymmetricSolidObject;

% dirichletBoundary2 = dirichletClass(dofObject);
% dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
% dirichletBoundary2.nodalDof = 2;
% dirichletBoundary2.timeFunction = @(t, XYZ) lengthY/2;
% dirichletBoundary2.masterObject = axisymmetricSolidObject;

% Nodal Forces
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = axisymmetricSolidObject;
nodalLoadObject.loadVector = M0/lengthX*[0; 1]*2*pi;
nodalLoadObject.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:,1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:,2) == lengthY/2);

nodalLoadObject2 = nodalLoadClass(dofObject);
nodalLoadObject2.masterObject = axisymmetricSolidObject;
nodalLoadObject2.loadVector = -M0/lengthX*[0; 1]*2*pi;
nodalLoadObject2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:,1) == innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:,2) == lengthY/2);

%% solver
dofObject = runNewton(setupObject, dofObject);
nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 1) - axisymmetricSolidObject.qN(nodeToMeasure, 1);
disp(['End displacement: ', num2str(endDisplacement)]);
exactDisplacement = uAnalytical;
disp(['Ratio FEM/Exact: ', num2str(endDisplacement/exactDisplacement)]);