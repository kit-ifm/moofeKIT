%% Dynamical Port-Hamiltonian Timoshenko beam with pointload at free end
clearvars; close all;

%% Parameters
lengthBeam = 1;
beamType = 'Timoshenko';
numberOfElements = 10;
E = 1840000;
R = 0.1;
nGP = 2;
rho = 920;
A = pi*R^2;
I = (pi*R^4)/4;
rhoA = rho*A;
shearCorFact = 6/7; %shear correction factor (rectangle 5/6, circle )

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-7;
setupObject.integrator = 'Midpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 1; %select 1 for bending moment and 2 for shear force
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
beamObject.elementDisplacementType = 'mixedPH';
beamObject.materialObject.E = E;
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
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
dirichletObject.nodalDof = [1, 2];
dirichletObject.masterObject = beamObject;

% neumann boundaries
Q = -rhoA*[5;0];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);
nodalLoadObject.timeFunction = str2func('@(t) sin(pi*t/0.2).*(t<=0.2+1e-12)');

%% solver
dofObject = runNewton(setupObject, dofObject);

% Get energy quantities
kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
time = zeros(setupObject.totalTimeSteps+1,1);

for j = 1:setupObject.totalTimeSteps+1
    time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
    kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
    potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.elementData(j).strainEnergy;
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