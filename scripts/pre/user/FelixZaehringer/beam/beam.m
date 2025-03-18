%% Aufgabe 13 a) - Grundlagen Finite Elemente

lengthBeam = 100;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'GFEAufgabe13a';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-6;
setupObject.integrator = 'Endpoint';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
beamObject = beamClass(dofObject);
beamObject.materialObject.name = 'Hooke';
beamObject.theory = 'timoshenko';
beamObject.elementDisplacementType = 'displacement';
beamObject.elementNameAdditionalSpecification = 'PetrovGalerkinANS';
widthBeam = 1 / 50 * lengthBeam;
beamObject.materialObject.E = 2.1 * 10^6;
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E / (2 * (1 + nu));
beamObject.materialObject.I = widthBeam^4 / 12;
beamObject.materialObject.A = widthBeam^2;
beamObject.materialObject.rho = 0;

numberOfElements = 3;

[nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
nodes = nodes + 1 / 2 * lengthBeam;
beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];

beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = 1;
beamObject.shapeFunctionObject.numberOfGausspoints = 2;

% dirichlet boundaries
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2];
dirichletObject.masterObject = beamObject;

% neumann boundaries
Q = 100;
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = -Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);

%% solver
dofObject = runNewton(setupObject, dofObject);
displacementLastNode = dofObject.listContinuumObjects{1}.qN1(end, 1);

analyticalSolution = -1/(beamObject.materialObject.E * beamObject.materialObject.I) * 1/3* Q * lengthBeam^3; % bernoulli beam theory

disp(['Displacement: ', num2str(displacementLastNode)])
disp(['Ratio: ', num2str(displacementLastNode/analyticalSolution)]);
