%% necessary input/ values
innerRadius = 4.12;
outerRadius = 4.32;
width = 0.1;
numberOfElements = 16;
nelZ=1;
E = 1e7;
nu = 0.3;
load = 'P1'; % P1 or P2

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
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Endpoint';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
order = 1;
serendipity = false;
solidObject.materialObject.name = 'Hooke';
%solidObject.elementDisplacementType = 'eas';
%solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
%solidObject.mixedFEObject.typeShapeFunctionData = 12;
solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinZhouEtAl2017'; % _Version_2_ für alte (verbesserte) Version , _Version_1_ für neu geschriebene Version , ohne alles für alte "schlechte" version
[solidObject.meshObject.nodes, solidObject.meshObject.edof, bounEdof] = meshCurvedBeam3D(innerRadius, outerRadius, numberOfElements, nelZ, width , order, serendipity);
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidObject.materialObject.mu = E/(2*(1+nu));
solidObject.materialObject.rho = 0;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = order;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:, 2) == 0);
dirichletBoundary.nodalDof = 1:3;
dirichletBoundary.masterObject = solidObject;

% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = solidObject;
if strcmpi(load, 'P1')
    neumannObject.loadVector = -[0; 1; 0]/0.1/0.2;
elseif strcmpi(load, 'P2')
    neumannObject.loadVector = -[0; 0; 1]/0.1/0.2;
end
neumannObject.meshObject.edof = bounEdof.SX1;

%% solver
dofObject = runNewton(setupObject, dofObject);
evalPointX = 0;

evalPointZ = 0.05;
evalPointY = outerRadius;
nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8 & abs(solidObject.meshObject.nodes(:, 3)-evalPointZ) < 1e-8);
if strcmpi(load, 'P1')
    endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
    normalizedEndDisplacement = round(endDisplacement/0.08734*(-1), 3);
elseif strcmpi(load, 'P2')
    endDisplacement = solidObject.qN1(nodeToMeasure, 3) - solidObject.qN(nodeToMeasure, 3);
    normalizedEndDisplacement = round(endDisplacement/0.5022*(-1), 3);
end
disp(['Normalized end displacement: ', num2str(normalizedEndDisplacement)]);