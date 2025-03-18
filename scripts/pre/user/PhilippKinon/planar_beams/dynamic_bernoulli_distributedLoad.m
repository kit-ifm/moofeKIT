clearvars; close all;

lengthBeam = 1;
beamType = 'Bernoulli';
numberOfElements = 10;
E = 1840000;
R = 0.1;
nGP = 2;
rho = 920;
A = pi*R^2;
I = (pi*R^4)/4;
rhoA = rho*A;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime = 10;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Midpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 2; %select 1 for bending moment and 2 for shear force
setupObject.plotObject.view = [0,90];
setupObject.plotObject.setXminXmax([-0.1;-0.7;-1],[1;0.7;1]);
setupObject.plotObject.border = 0.1;
setupObject.plotObject.plotInitialConfig = true;
setupObject.plotObject.showGrid = true;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
beamObject = beamClass(dofObject);
beamObject.materialObject.name = "Hooke";
beamObject.theory = beamType;
beamObject.elementDisplacementType = 'displacement';
beamObject.materialObject.E = E;
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
beamObject.materialObject.I = I;
beamObject.materialObject.A = A;
beamObject.materialObject.rho = rho;

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

% distributed load
q = [rhoA*9.81;0];
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
