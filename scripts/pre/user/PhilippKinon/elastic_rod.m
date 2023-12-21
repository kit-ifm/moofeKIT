%% Elastic rod
%
% based on: oneDimensionalContinuum.m
%
% problem: hyperelastic or linear elastic rod, fixed on one end, 
%          time-dependent force on the other
% spatial discretization: linear Lagrange shapefunction for disp
%                         discont. constant ansatz for strain, 
%                         mixed formulation
% time discretization: midpoint rule or discrete gradient
%
% user: PK

%% Setup
setupObject = setupClass;
setupObject.saveObject.fileName = 'rod';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 50;
setupObject.totalTime = 5;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = 1;
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Midpoint'; 
%setupObject.integrator = 'DiscreteGradient';
dofObject = dofClass;   % required object for dof and object handling

%% Continuum object
solidObject = solidClass(dofObject);
solidObject.materialObject.name = 'Hooke';
%solidObject.materialObject.name = 'Hyperelastic';
%solidObject.elementDisplacementType = 'displacement';
solidObject.elementDisplacementType = 'mixedPH'; %%pH Formulation with Green-Lagrangian strains
solidObject.elementNameAdditionalSpecification = 'E';

E = 300000;
nu = 3e-01;
solidObject.materialObject.E = E;
solidObject.materialObject.nu = nu;
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));                 % Lame parameter
solidObject.materialObject.mu = E/(2*(1+nu));
solidObject.materialObject.rho = 1500; 
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = false;

%% Geometry
numberOfElementsOfCrime = 4;
length = 2;
order  = 1;
number_of_gausspoints = 2;

% Mesh
[solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = meshOneDimensional(length,numberOfElementsOfCrime,order);
solidObject.meshObject.nodes = solidObject.meshObject.nodes + length/2; %% TO DO: write adapted mesh

% Shapefunctions
solidObject.dimension = 1;
solidObject.shapeFunctionObject.order = order;
solidObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;

% Dirichlet BC (fixed at x=0)
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1;
dirichletObject.masterObject = solidObject; %% TO Do: use a term which is not racist
dirichletObject.timeFunction = str2func('@(t,X) 0*X');

% neumann boundary conditions
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = solidObject;
nodalLoadObject.loadVector = 10000;
nodalLoadObject.timeFunction = str2func('@(t)  sin(pi/2 * t).*(t<=1+1e-12)');
nodalLoadObject.nodeList= size(solidObject.meshObject.nodes,1);
loadingTime = 1; % Time until loading happens

%% Solver
dofObject = runNewton(setupObject,dofObject);

%% Postprocessing
% Plot the problem
plot(solidObject,setupObject)

% Get energy quantities
kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
time = zeros(setupObject.totalTimeSteps+1,1);

for j = 1:setupObject.totalTimeSteps
    time(j+1) = j*setupObject.totalTime/setupObject.totalTimeSteps;
    kineticEnergy(j+1) = dofObject.postDataObject.energyJournal(j+1).EKin;
    potentialEnergy(j+1) = dofObject.listContinuumObjects{1,1}.ePot(j+1).strainEnergy;
end

% Compute increment
figure()
plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
legend('T', 'V', 'H');

% Plot increment
figure()
energyIncrement = zeros(setupObject.totalTimeSteps,1);
timeWithoutLoading = loadingTime:setupObject.totalTime/setupObject.totalTimeSteps:setupObject.totalTime;

for i = 1:setupObject.totalTimeSteps
    energyIncrement(i) = (kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)); 
end

plot(timeWithoutLoading(1:end-1),energyIncrement(loadingTime/setupObject.totalTime*setupObject.totalTimeSteps+1:end))
legend('Energiedifferenz');
