%% Dynamical GE beam flying spaghetti
clearvars; close all;

%% Parameters
lengthBeam = 10;
order  = 1;
beamType = 'GeometricallyExact';
numberOfElements = 10;
EA = 10000;
kGA = 10000;
nGP = 2;
rhoA = 1;
EI = 100;
rhoI = 10;
rho = 1; %dummy
stepsize=0.1;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTime = 15;
setupObject.totalTimeSteps = setupObject.totalTime/stepsize;
setupObject.newton.tolerance = 1e-11;
setupObject.integrator = 'Midpoint';
setupObject.plotObject.flag = true;
setupObject.plotObject.steps = 5;
setupObject.plotObject.keepFormerPlots = true;
setupObject.plotObject.postPlotType = "stress";
setupObject.plotObject.stress.component = 3; %select 3 for bending moment, 1 for normal force, and 2 for shear force
setupObject.plotObject.view = [0,90];
setupObject.plotObject.lineWeight = 1.1;
setupObject.plotObject.setXminXmax([-1;-1;-1],[18;9;1]);
setupObject.plotObject.border = 0;
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
beamObject.materialObject.EA = EA;
beamObject.materialObject.kGA = kGA;
beamObject.materialObject.EI = EI;
beamObject.materialObject.rhoA = rhoA;
beamObject.materialObject.rhoI = rhoI;
beamObject.materialObject.rho = rho;

%% Spatial discretiaztion
startpoint = [0;8];
endpoint = [6;0];
length = lengthBeam;

% Mesh
[beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements,order,startpoint,endpoint);
beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = order;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;

% neumann boundaries
P = 8;
T = 80;
F = [P;0;-T];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = F;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == endpoint(1));
nodalLoadObject.timeFunction = str2func('@(t) (t<=2.5+1e-12)');

%% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;
dofObject = runNewton(setupObject, dofObject);

% matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'snapshots.tikz', 'showInfo', false, 'floatformat', '%.4g')


% Get energy quantities
kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
time = zeros(setupObject.totalTimeSteps+1,1);

for j = 1:setupObject.totalTimeSteps+1
    time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
    kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
    potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).strainEnergy;
end

figure()
plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
legend('T', 'V', 'H');
%matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'energy.tikz', 'showInfo', false, 'floatformat', '%.4g')


% Angular momentum
figure()
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',1);
plot(time,angularMomentum)
legend('Total Angular Momentum');
%matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'ang_mom.tikz', 'showInfo', false, 'floatformat', '%.4g')


% Plot increment
figure()
energyIncrement = zeros(setupObject.totalTimeSteps,1);
angMomIncrement = zeros(setupObject.totalTimeSteps,1);

for i = 1:setupObject.totalTimeSteps
    energyIncrement(i) = abs((kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)));
    angMomIncrement(i) = abs(angularMomentum(i+1) - angularMomentum(i));
end
plot(time(1:end-1),energyIncrement)
legend('Energy increment');
%matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'energy_increment.tikz', 'showInfo', false, 'floatformat', '%.4g')

figure()
plot(time(1:end-1),angMomIncrement)
legend('Angular momentum increment');
%matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'ang_mom_increment.tikz', 'showInfo', false, 'floatformat', '%.4g')
