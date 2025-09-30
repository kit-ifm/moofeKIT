%% Dynamical GE beam test
clearvars;

%% Parameters
order  = 1;
beamType = 'GeometricallyExact';
nGP = 2;

lengthBeam1 = 0.5;
area = 7.854e-5;
I = 4.909e-10;
E2 = 0.5e7;
rho = 2770;
nu = 0.3;
G2 = E2/(2*(1+nu));

numberOfElements = 4;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
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
beamObject.materialObject.E = E2;
beamObject.materialObject.A = area;
beamObject.materialObject.I = I;
beamObject.materialObject.G = G2;
beamObject.materialObject.rho = rho;
beamObject.materialObject.shearCorrectionCoefficient = 1;

% Spatial discretiaztion
startpoint = [0;0];
endpoint = [lengthBeam1;0];
length = lengthBeam1;

% Mesh
[beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements,order,startpoint,endpoint);
beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = order;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
numberOfDispDOFFirstBeam = 3 * max(max(beamObject.meshObject.edof)); % displacements dofs
numberOfDOFFirstBeam = numberOfDispDOFFirstBeam + numberOfElements * 3; % strain dofs

% dirichlet boundaries
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(beamObject.meshObject.nodes(:, 1) == startpoint(1));
dirichletObject1.nodalDof = [1, 2];
dirichletObject1.masterObject = beamObject;

%Bodyforce
bodyForceObject = bodyForceClass(dofObject);
bodyForceObject.typeOfLoad = 'deadLoad';
bodyForceObject.masterObject = beamObject;
bodyForceObject.loadFunction = rho*area*[0; -9.81];
bodyForceObject.dimension = 2;
bodyForceObject.timeFunction = str2func('@(t) 1');
bodyForceObject.meshObject.edof = beamObject.meshObject.edof;
bodyForceObject.shapeFunctionObject.order = 1;
bodyForceObject.shapeFunctionObject.numberOfGausspoints = nGP;

% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;


%% Solve!
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

figure()
plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
legend('T', 'V', 'H');

% Angular momentum
figure()
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',1);
plot(time,angularMomentum)
legend('Total Angular Momentum');

% Plot increment
figure()
energyIncrement = zeros(setupObject.totalTimeSteps,1);
angMomIncrement = zeros(setupObject.totalTimeSteps,1);

for i = 1:setupObject.totalTimeSteps
    energyIncrement(i) = abs((kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)));
    angMomIncrement(i) = angularMomentum(i+1) - angularMomentum(i);
end
plot(time(1:end-1),energyIncrement)
legend('Energy increment');

figure()
plot(time(1:end-1),angMomIncrement)
legend('Angular momentum increment');
