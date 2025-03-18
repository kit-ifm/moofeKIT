%% Static GE beam with point torque at free end
clearvars; close all;

%% Parameters
lengthBeam = 1;
beamType = 'GeometricallyExact';
numberOfElements = 8;
E = 184000;
R = 0.186;
nGP = 2;
rho = 0;
nu = 0.3;
G = E/(2*(1+nu));
A = pi*R^2;
I = (pi*R^4)/4;
shearCorFact = 1; %shear correction factor (rectangle 5/6, circle 6/7)

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1000;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-11;
setupObject.integrator = 'Endpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 3; %select 1 for bending moment and 2 for shear force
setupObject.plotObject.view = [0,90];
setupObject.plotObject.setXminXmax([-0.3;-0.1;-1],[1;0.7;1]);
setupObject.plotObject.border = 0.1;
setupObject.plotObject.plotInitialConfig = true;
setupObject.plotObject.showGrid = true;
setupObject.plotObject.colorscheme = 'viridis';

dofObject = dofClass; % required object for dof and object handling
dofObject.postDataObject.storeStateFlag = true;

%% continuum Objects
beamObject = beamClass(dofObject,2);
beamObject.materialObject.name = "Hooke";
beamObject.theory = beamType;
beamObject.elementDisplacementType = 'mixedPH';
beamObject.materialObject.E = E;
beamObject.materialObject.G = G;
beamObject.materialObject.I = I;
beamObject.materialObject.A = A;
beamObject.materialObject.rho = rho;
beamObject.materialObject.shearCorrectionCoefficient = shearCorFact;

[nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
nodes = nodes + 1 / 2 * lengthBeam;
beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];

beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = 1;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;

% dirichlet boundaries
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2, 3];
dirichletObject.masterObject = beamObject;

% neumann boundaries
Q = [0;0;2*E*I*pi/lengthBeam];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);
nodalLoadObject.timeFunction = str2func('@(t) t');

%% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;
dofObject = runNewton(setupObject, dofObject);