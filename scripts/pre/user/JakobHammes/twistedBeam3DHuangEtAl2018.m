%% necessary input/ values
beamlength = 12;
height = 1.1;
width = 0.32; % Zhou
%width = 0.0032; % Huang
angle = 90; % degrees
nelX = 4;
nelY = 2;
nelZ = 1;
E = 2.9e7;
nu = 0.22;

loading = 'P1'; % in plane
%loading = 'P2'; % out of plane


%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'twistedBeam3D';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 11;
setupObject.plotObject.view = 2;
setupObject.usePreconditioning = false;
setupObject.newton.tolerance = 1e-2;
setupObject.integrator = 'Endpoint';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidShellObject = solidShellClass(dofObject);
order = 1;
serendipity = false;
solidShellObject.materialObject.name = 'Hooke';
solidShellObject.elementNameAdditionalSpecification = 'PetrovGalerkinHuangEtAl2018';
[solidShellObject.meshObject.nodes, solidShellObject.meshObject.edof, bounEdof] = meshTwistedBeam3D(beamlength, height, width, angle, nelX, nelY, nelZ, order, serendipity);solidShellObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidShellObject.materialObject.mu = E/(2*(1+nu));
solidShellObject.materialObject.rho = 0;
solidShellObject.dimension = 3;
solidShellObject.shapeFunctionObject.order = order;
solidShellObject.shapeFunctionObject.numberOfGausspoints = 8;

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(solidShellObject.meshObject.nodes(:, 1) == 0);
dirichletBoundary.nodalDof = 1:3;
dirichletBoundary.masterObject = solidShellObject;

% Neumann boundary
F = 1;
%F = 1/10e6; % Belytschko
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = solidShellObject;

switch loading
    case 'P1' % in plane
        neumannObject.loadVector = [0; 0; F]/width/height;
    case 'P2' % out of plane
        neumannObject.loadVector = -[0; F; 0]/width/height;
end
neumannObject.meshObject.edof = bounEdof.SX2;

%% solver
dofObject = runNewton(setupObject, dofObject);
evalPointX = beamlength;
evalPointZ = 0;
evalPointY = width/2;
nodeToMeasure = find(abs(solidShellObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(solidShellObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8 & abs(solidShellObject.meshObject.nodes(:, 3)-evalPointZ) < 1e-8);

switch loading
    case 'P1' % in plane
        endDisplacement = solidShellObject.qN1(nodeToMeasure, 3) - solidShellObject.qN(nodeToMeasure, 3);
        standardValue = 0.005424; % from Zhou 2017 for width = 0.32
        %standardValue = 0.005256; % Huang?/Belytschko for width = 0.0032
    case 'P2' % out of plane
        endDisplacement = solidShellObject.qN1(nodeToMeasure, 2) - solidShellObject.qN(nodeToMeasure, 2);
        standardValue = 0.001754; % from Zhou 2017 for width = 0.32
        %standardValue = 0.001294; % Huang?/Belytschko for width = 0.0032
end
normalizedEndDisplacement = endDisplacement/standardValue; 

disp(['Max. displacement: ', num2str(endDisplacement)]);
disp(['Normalized end displacement: ', num2str(round(normalizedEndDisplacement, 4))]);