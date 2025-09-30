%% necessary input/ values
radius = 10;
thickness = 1;
innerRadius = radius - thickness/2;
outerRadius = radius + thickness/2;
width = 1;
numberOfElements = 10;
nelZ=1;
E = 1e7;
nu = 0.3;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'curvedBeam3D';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = 2;
setupObject.usePreconditioning = false;
setupObject.newton.tolerance = 1e-3;
setupObject.integrator = 'Endpoint';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidShellObject = solidShellClass(dofObject);
order = 1;
serendipity = false;
solidShellObject.materialObject.name = 'Hooke';
solidShellObject.elementNameAdditionalSpecification = 'PetrovGalerkinHuangEtAl2018';
[solidShellObject.meshObject.nodes, solidShellObject.meshObject.edof, bounEdof] = meshCurvedBeam3D(innerRadius, outerRadius, numberOfElements, nelZ, width , order, serendipity);
solidShellObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidShellObject.materialObject.mu = E/(2*(1+nu));
solidShellObject.materialObject.rho = 0;
solidShellObject.dimension = 3;
solidShellObject.shapeFunctionObject.order = order;
solidShellObject.shapeFunctionObject.numberOfGausspoints = 8;

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(solidShellObject.meshObject.nodes(:, 2) == 0);
dirichletBoundary.nodalDof = 1:3;
dirichletBoundary.masterObject = solidShellObject;

% Neumann boundary
F = 1;
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = solidShellObject;
neumannObject.loadVector = -[0; F; 0]/width/thickness;
neumannObject.meshObject.edof = bounEdof.SX1;

%% solver
dofObject = runNewton(setupObject, dofObject);
evalPointX = 0;
evalPointZ = -width/2;
evalPointY = innerRadius;
nodeToMeasure = find(abs(solidShellObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(solidShellObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8 & abs(solidShellObject.meshObject.nodes(:, 3)-evalPointZ) < 1e-8);
endDisplacement = solidShellObject.qN1(nodeToMeasure, 2) - solidShellObject.qN(nodeToMeasure, 2);
analyticalSolution = 3*pi*F*radius^3/(E*width*thickness^3);
normalizedEndDisplacement = endDisplacement/analyticalSolution*(-1);
disp(['Normalized end displacement: ', num2str(round(normalizedEndDisplacement, 4))]);