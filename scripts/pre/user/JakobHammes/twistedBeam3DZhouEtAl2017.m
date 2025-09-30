%% necessary input/ values
beamlength = 12;
height = 1.1;
width = 0.32;   % Zhou
%width = 0.0032;  % Huang
angle = 90; % degrees
nelX = 12;
nelY = 2;
nelZ = 1;
order = 1;
serendipity = false;
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
setupObject.plotObject.flag = false;
setupObject.plotObject.stress.component = 11;
setupObject.plotObject.view = 2;
setupObject.plotObject.postPlotType = 'stress';
%setupObject.plotObject.stress.component = 'maxPrincipalStress'; 
setupObject.newton.tolerance = 1e-2;
setupObject.integrator = 'Endpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.materialObject.name = 'Hooke';
solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinZhouEtAl2017';
%solidObject.elementDisplacementType = 'eas';
%solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';

solidObject.mixedFEObject.typeShapeFunctionData = 12;
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidObject.materialObject.mu = E/(2*(1+nu));
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = E;
solidObject.materialObject.nu = nu;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = order;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;
[solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEdof] = meshTwistedBeam3D(beamlength, height, width, angle, nelX, nelY, nelZ, order, serendipity);

% Dirichlet boundary
dirichletBoundary = dirichletClass(dofObject);
dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletBoundary.nodalDof = 1:3;
dirichletBoundary.masterObject = solidObject;

% Neumann boundary
F = 1;
%F = 1*10e-6; % Belytschko
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'area';
neumannObject.masterObject = solidObject;

switch loading
    case 'P1' % in plane
        neumannObject.loadVector = [0; 0; F]/width/height;
    case 'P2' % out of plane
        neumannObject.loadVector = -[0; F; 0]/width/height;
end
neumannObject.meshObject.edof = bounEdof.SX2;


%% solver
dofObject = runNewton(setupObject,dofObject);
nodeToMeasure = find(solidObject.meshObject.nodes(:, 1) == beamlength & solidObject.meshObject.nodes(:, 2) == width/2 & solidObject.meshObject.nodes(:, 3) == 0);

switch loading
    case 'P1' % in plane
        endDisplacement = solidObject.qN1(nodeToMeasure, 3) - solidObject.qN(nodeToMeasure, 3);
        standardValue = 0.005424; % from Zhou 2017 for width = 0.32
        %standardValue = 0.005256; % Huang?/Belytschko for width = 0.0032
    case 'P2' % out of plane
        endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
        standardValue = 0.001754; % from Zhou 2017 for width = 0.32
        %standardValue = 0.001294; % Huang?/Belytschko for width = 0.0032
end
normalizedEndDisplacement = endDisplacement/standardValue; 

disp(['Max. displacement: ', num2str(endDisplacement)]);
disp(['Normalized end displacement: ', num2str(round(normalizedEndDisplacement, 4))]);
