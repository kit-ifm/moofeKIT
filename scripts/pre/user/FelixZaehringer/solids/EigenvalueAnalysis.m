%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'block';
setupObject.saveObject.saveData = false;
setupObject.totalTime = 2;
setupObject.plotObject.flag = true;
setupObject.newton.tolerance = 1e-5;
setupObject.integrator = 'Endpoint';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.dimension = 2;
solidObject.elementDisplacementType = 'displacement';
% solidObject.elementDisplacementType = 'pianSumihara';
% solidObject.elementDisplacementType = 'selectiveReducedIntegration';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
solidObject.shapeFunctionObject.order = 1;
serendipity = false;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;
solidObject.materialObject.name = 'HookeEVZ';
% solidObject.materialObject.name = 'Hooke2D';
% solidObject.materialObject.name = 'NeoHooke2D';
solidObject.mixedFEObject.condensation = true;
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1e4;
solidObject.materialObject.nu = 0.2;
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(2, 2, 1, 1, solidObject.shapeFunctionObject.order, serendipity);
% solidObject.meshObject.nodes(:, :) = solidObject.meshObject.nodes(:, :) + 25;

% dirchlet boundary conditions
% dirichletObject1 = dirichletClass(dofObject);
% dirichletObject1.masterObject = solidObject;
% dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
% dirichletObject1.nodalDof = 1;
% 
% dirichletObject = dirichletClass(dofObject);
% dirichletObject.masterObject = solidObject;
% dirichletObject.nodeList = find((solidObject.meshObject.nodes(:, 1) == 0) & (abs(solidObject.meshObject.nodes(:, 2)) == 0));
% dirichletObject.nodalDof = 1:2;

% neumann boundary conditions
% neumannObject = neumannClass(dofObject);
% neumannObject.loadGeometry = 'line';
% neumannObject.masterObject = solidObject;
% neumannObject.loadVector = [-200; 0];
% neumannObject.timeFunction = @(t) sin(pi*t/2) .* (t <= 1);
% neumannObject.meshObject.edof = edofBoundary.SX2;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
% yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
% fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)

%% postprocessing - eigenvalue analysis
K = full(dofObject.K);
% obtain eigenvalues/vectors
[eigVectors, eigValues] = eig(K);
eigValues = eigValues*ones(size(eigValues, 1), 1);

% sort eigenvalues/vectors
[eigValues, I] = sort(eigValues);
eigVectors = eigVectors(:, I);

% plot eigenvectors
% figure;
% plot(solidObject, setupObject);

% disp eigenvalues
disp(eigValues);
