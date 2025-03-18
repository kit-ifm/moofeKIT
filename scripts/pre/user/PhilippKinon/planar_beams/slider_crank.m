%% Dynamical GE beam multibody system
clearvars; close all;

%% Parameters
order  = 1;
beamType = 'GeometricallyExact';
nGP = 2;

% data from stander and stein
lengthBeam1 = 0.1524;
lengthBeam2 = 0.3048;
diameter = 0.00635;
area = pi*diameter^2/4;
I = (pi*(diameter/2)^4)/4;
E1 = 20684271840000;
E2 = 206842718400;
rho = 7834.6397;
nu = 0.0;
G1 = E1/(2*(1+nu));
G2 = E2/(2*(1+nu));

numberOfElements1 = 3;
numberOfElements2 = 8;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 400;
setupObject.totalTime = 10;
setupObject.newton.tolerance = 1e-10;
setupObject.integrator = 'Midpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 3; %select 3 for bending moment, 1 for normal force, and 2 for shear force
setupObject.plotObject.view = [0,90];
setupObject.plotObject.setXminXmax([-0.2;-0.35;-1],[0.5;0.35;1]);
setupObject.plotObject.border = 0.1;
setupObject.plotObject.plotInitialConfig = true;
setupObject.plotObject.showGrid = true;
setupObject.plotObject.colorscheme = 'viridis';

dofObject = dofClass; % required object for dof and object handling
dofObject.postDataObject.storeStateFlag = true;

%% beam 1
beamObject = beamClass(dofObject,2);
beamObject.materialObject.name = "Hooke";
beamObject.theory = beamType;
beamObject.elementDisplacementType = 'mixedPH';
beamObject.materialObject.E = E1;
beamObject.materialObject.A = area;
beamObject.materialObject.I = I;
beamObject.materialObject.G = G1;
beamObject.materialObject.rho = rho;
beamObject.materialObject.shearCorrectionCoefficient = 1;

% Spatial discretiaztion
startpoint = [0;0];
endpoint = [lengthBeam1;0];
length = lengthBeam1;

% Mesh
[beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements1,order,startpoint,endpoint);
beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = order;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
numberOfDispDOFFirstBeam = 3 * max(max(beamObject.meshObject.edof)); % displacements dofs
numberOfDOFFirstBeam = numberOfDispDOFFirstBeam + numberOfElements1 * 3; % strain dofs

% dirichlet boundaries
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(beamObject.meshObject.nodes(:, 1) == startpoint(1));
dirichletObject1.nodalDof = [1, 2];
dirichletObject1.masterObject = beamObject;

% neumann boundaries
M = 0.01;
F = [0;0;M];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = F;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == startpoint(1));
nodalLoadObject.timeFunction = str2func('@(t) (t<=1+1e-12)');

% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;

%% beam 2
beamObject2 = beamClass(dofObject,2);
beamObject2.materialObject.name = "Hooke";
beamObject2.theory = beamType;
beamObject2.elementDisplacementType = 'mixedPH';
beamObject2.materialObject.E = E2;
beamObject2.materialObject.A = area;
beamObject2.materialObject.I = I;
beamObject2.materialObject.G = G2;
beamObject2.materialObject.rho = rho;
beamObject2.materialObject.shearCorrectionCoefficient = 1;

% Spatial discretiaztion
startpoint = [lengthBeam1;0];
endpoint = [lengthBeam1+lengthBeam2;0];
length = lengthBeam2;

% Mesh
[beamObject2.meshObject.nodes,beamObject2.meshObject.edof,edofNeumann2] = linearString(length,numberOfElements2,order,startpoint,endpoint);
beamObject2.dimension = 1;
beamObject2.shapeFunctionObject.order = order;
beamObject2.shapeFunctionObject.numberOfGausspoints = nGP;

% dirichlet boundaries
dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = find(beamObject2.meshObject.nodes(:, 1) == endpoint(1));
dirichletObject2.nodalDof = 2;
dirichletObject2.masterObject = beamObject2;

%% solver
beamObject2.numericalTangentObject.computeNumericalTangent = true;
beamObject2.numericalTangentObject.showDifferences = false;

%% Include constraints
constraintObject1 = constraintClass(dofObject);
constraintObject1.constrainedElement = numberOfElements1 + 1;
constraintObject1.constrainedDOFinElement = [numberOfDOFFirstBeam + 1, numberOfDOFFirstBeam + 2];
constraintObject1.constrainedDOFinOtherElement = [numberOfDispDOFFirstBeam - 2 , numberOfDispDOFFirstBeam - 1];

%% Solve!
dofObject = runNewton(setupObject, dofObject);

% Get energy quantities
kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
time = zeros(setupObject.totalTimeSteps+1,1);

for j = 1:setupObject.totalTimeSteps+1
    time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
    kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
    potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).strainEnergy + dofObject.listContinuumObjects{1,4}.ePot(j).strainEnergy;
end

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
