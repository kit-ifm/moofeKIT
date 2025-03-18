%% GEB with pointload at free end (rotated to test objectivity)
clearvars;
%close all;

rotationAngle = pi/6;

%% Parameters
lengthBeam = 1;
beamType = 'GeometricallyExact';
numberOfElements = 10;
E = 184000;
R = 0.186;
nGP = 2;
rho = 0;
A = pi*R^2;
I = (pi*R^4)/4;
rhoA = rho*A;
shearCorFact = 6/7; %shear correction factor (rectangle 5/6, circle 6/7)

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-10;
setupObject.integrator = 'Midpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 3; %select 3 for bending moment, 1 for normal force, and 2 for shear force
setupObject.plotObject.view = [0,90];
setupObject.plotObject.setXminXmax([-0.1;-0.7;-1],[1;0.7;1]);
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
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
beamObject.materialObject.I = I;
beamObject.materialObject.A = A;
beamObject.materialObject.rho = rho;
beamObject.materialObject.shearCorrectionCoefficient = shearCorFact;

startpoint = [0;0];
endpoint = [cos(rotationAngle)*lengthBeam;sin(rotationAngle)*lengthBeam];
length = lengthBeam;

% Mesh
[beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements,1,startpoint,endpoint);
beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = 1;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;

% dirichlet boundaries
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2, 3];
dirichletObject.masterObject = beamObject;

% neumann boundaries
Q = 1000*[-sin(rotationAngle);cos(rotationAngle);0];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == endpoint(1));
nodalLoadObject.timeFunction = str2func('@(t) t');

%% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;
dofObject = runNewton(setupObject, dofObject);

% Get energy quantities
kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
time = zeros(setupObject.totalTimeSteps+1,1);

for j = 1:setupObject.totalTimeSteps+1
    time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
    kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
    potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).strainEnergy;
end

% Compute increment
figure()
plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
legend('T', 'V', 'H');

% Plot increment
figure()
energyIncrement = zeros(setupObject.totalTimeSteps,1);

for i = 1:setupObject.totalTimeSteps
    energyIncrement(i) = abs((kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)));
end
plot(time(1:end-1),energyIncrement)
legend('Energy increment');