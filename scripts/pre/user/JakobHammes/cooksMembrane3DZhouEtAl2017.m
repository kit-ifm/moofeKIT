%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembrane3D';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = 2;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = 'maxPrincipalStress'; 
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Endpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.materialObject.name = 'Hooke';
solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinZhouEtAl2017';
% solidObject.elementDisplacementType = 'eas';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
% solidObject.mixedFEObject.typeShapeFunctionData = 12;

E = 1;
nu = 1/3;
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidObject.materialObject.mu = E/(2*(1+nu));
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = E;
solidObject.materialObject.nu = nu;
nel = 2;
nelZ = 1;
[solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,nelZ,48,44,16,16);
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;

% Dirichlet boundary
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1:2;
dirichletObject.masterObject = solidObject;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0 & solidObject.meshObject.nodes(:,2)==0 & solidObject.meshObject.nodes(:,3)==0);
dirichletObject.nodalDof = 3;
dirichletObject.masterObject = solidObject;

neumannObject = neumannClass(dofObject);
neumannObject.masterObject = solidObject;
neumannObject.loadGeometry = 'area';
neumannObject.loadVector = [0;1;0]/16;
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
nodeToMeasure = find(solidObject.meshObject.nodes(:, 1) == 48 & solidObject.meshObject.nodes(:, 2) == 52 & solidObject.meshObject.nodes(:, 3) == 0);
disp(['Max. displacement: ', num2str(solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2))]);

% stress analysis
if strcmpi(setupObject.plotObject.stress.component, 'maxPrincipalStress')
    nodeToMeasure = find(abs(solidObject.meshObject.nodes(:,1) - 24) < 1e-8 & abs(solidObject.meshObject.nodes(:,2) - 22) < 1e-8 & solidObject.meshObject.nodes(:, 3) == 0);
    stress = dofObject.listContinuumObjects{1}.plotData.nodalColorData(nodeToMeasure);
    disp(['Max. principle stress (Point A): ', num2str(stress)]);
elseif strcmpi(setupObject.plotObject.stress.component, 'minPrincipalStress')
    nodeToMeasure = find(abs(solidObject.meshObject.nodes(:,1) - 24) < 1e-8 & abs(solidObject.meshObject.nodes(:,2) - 52) < 1e-8 & solidObject.meshObject.nodes(:, 3) == 0);
    stress = dofObject.listContinuumObjects{1}.plotData.nodalColorData(nodeToMeasure);
    disp(['Min. principle stress (Point B): ', num2str(stress)]);
end
