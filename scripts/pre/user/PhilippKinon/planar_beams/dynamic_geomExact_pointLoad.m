%% GEB with pointload at free end
clearvars;
%close all;

rotationAngle = 0;
class = 'beam'; % beam, beamVelocity

%% Parameters
lengthBeam = 1;
beamType = 'GeometricallyExact';
numberOfElements = 10;
E = 184000;
R = 0.186;
nGP = 2;
selRednGP = 1;
rho = 920;
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
setupObject.plotObject.postPlotType = "time";
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
if strcmpi(class, 'beam')
    beamObject = beamClass(dofObject,2);
elseif strcmpi(class, 'beamVelocity')
    solidObject = beamVelocityClass(dofObject,2);
end


beamObject.materialObject.name = "Hooke";
beamObject.theory = beamType;
beamObject.elementDisplacementType = 'mixedPH';
%beamObject.elementDisplacementType = 'displacement';
%beamObject.elementNameAdditionalSpecification = 'StanderStein';
%beamObject.elementNameAdditionalSpecification = 'PhIrreducible';
beamObject.materialObject.E = E;
nu = 0.3;
beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
beamObject.materialObject.I = I;
beamObject.materialObject.A = A;
beamObject.materialObject.rho = rho;
beamObject.materialObject.shearCorrectionCoefficient = shearCorFact;

if strcmpi(class, 'beam')
    beamObject.materialObject.rho = rho;
elseif strcmpi(class, 'beamVelocity')
    beamObject.materialObject.rho = 0;
    solidObject.materialObject.rhoNew = rho;
end

startpoint = [0;0];
endpoint = [cos(rotationAngle)*lengthBeam;sin(rotationAngle)*lengthBeam];
length = lengthBeam;

% Mesh
[beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements,1,startpoint,endpoint);
beamObject.dimension = 1;
beamObject.shapeFunctionObject.order = 1;
beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
beamObject.selectiveReducedShapeFunctionObject.order = 1;
beamObject.selectiveReducedShapeFunctionObject.numberOfGausspoints = selRednGP;

% dirichlet boundaries
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2, 3];
dirichletObject.masterObject = beamObject;

if strcmpi(class, 'beamVelocity')
    %% TO CHECK: do sth else!?
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.masterObject = beamObject;
    dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = [1, 2, 3];
end

% neumann boundaries
Q = 1000*[-sin(rotationAngle);cos(rotationAngle);0];
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = beamObject;
nodalLoadObject.loadVector = Q;
nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == endpoint(1));
nodalLoadObject.timeFunction = str2func('@(t) sin(pi*t/0.2).*(t<=0.2+1e-12)');

%% solver
beamObject.numericalTangentObject.computeNumericalTangent = true;
beamObject.numericalTangentObject.showDifferences = false;
dofObject = runNewton(setupObject, dofObject);

% Get postprocessing quantities
time = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
potentialEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',1);
drift = getElementData(dofObject.postDataObject,dofObject,setupObject,'drift');

% Compute increments
energyIncrement = zeros(setupObject.totalTimeSteps,1);
angMomIncrement = zeros(setupObject.totalTimeSteps,1);
for i = 1:setupObject.totalTimeSteps
    energyIncrement(i) = abs((kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)));
    angMomIncrement(i) = abs(angularMomentum(i+1) - angularMomentum(i));
end

%Plot energies
figure()
plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
legend('T', 'V', 'H');

% Plot energy increment
figure()
plot(time(1:end-1),energyIncrement)
legend('Energy increment');

% Plot angular momentum
figure()
plot(time,angularMomentum)
legend('Total Angular Momentum');

% Plot angular momentum increment
figure()
plot(time(1:end-1),angMomIncrement)
legend('Angular momentum increment');

% Plot drift
figure()
plot(time,drift)
legend('Difference of independent strains and displacement-based computation');