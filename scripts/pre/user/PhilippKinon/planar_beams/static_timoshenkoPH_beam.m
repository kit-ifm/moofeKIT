clearvars;
%% Aufgabe 13 a) - Grundlagen Finite Elemente

lengthBeam = 100;
widthToLength = 1 / 50;
beamType = 'Timoshenko';
numberOfElements = 3;
nGP = 2;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-8;
setupObject.integrator = 'Endpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 1; %select 1 for bending moment and 2 for shear force
setupObject.plotObject.view = [0,90];

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
beamObject = beamClass(dofObject);
beamObject.materialObject.name = "Hooke";
beamObject.theory = beamType;
beamObject.elementDisplacementType = 'mixedPH'; 
widthBeam = widthToLength * lengthBeam;
beamObject.materialObject.E = 2.1 * 10^6;
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
beamObject.materialObject.I = widthBeam^4 / 12;
beamObject.materialObject.A = widthBeam^2;
beamObject.materialObject.rho = 0;
beamObject.materialObject.shearCorrectionCoefficient = 5/6;

[nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
nodes = nodes + 1 / 2 * lengthBeam;
beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];

beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = 1;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;

% dirichlet boundaries
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2];
dirichletObject.masterObject = beamObject;

% neumann boundaries
Q = [-100;0];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);

% distributed load
q = [-1;-10];
distributedForceObject = bodyForceClass(dofObject);
distributedForceObject.typeOfLoad = 'deadLoad';
distributedForceObject.masterObject = beamObject;
distributedForceObject.loadFunction = q;
distributedForceObject.timeFunction = str2func('@(t) 1');
distributedForceObject.dimension = 2; % second dimension required to fix bodyforce vector
distributedForceObject.meshObject.edof = beamObject.meshObject.edof;
distributedForceObject.shapeFunctionObject.order = 1;
distributedForceObject.shapeFunctionObject.numberOfGausspoints = nGP;

%% solver
dofObject = runNewton(setupObject, dofObject);
displacementLastNode = dofObject.listContinuumObjects{1}.qN1(end, 1);
