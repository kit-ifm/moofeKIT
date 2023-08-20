%% Elastic string
%
% based on: elastic_string.m
%
% problem: hyperelastic or linear elastic string, fixed on one end, time-dependent force on the other
% spatial discretization: linear Lagrange shapefunction, disp.-based
%                         formulation
% time discretization: midpoint rule or discrete gradient
%
% user: Philipp Kinon
% date: May 02, 2023
clearvars; close all;

%% Main parameters
totalTime = 1;
timeStepSize = 1e-2;
newtonTolerance = 1e-7;
integrator = 'Midpoint';
%integrator = 'DiscreteGradient';
%material = 'StVenant';
material = 'Hyperelastic';
EA = 10;
nu = 3e-01;
rhoA = 1;
numberOfElementsOfCrime = 30;
postprocess = true;
exportFolder = './';
nDIM = 2;

%% Simulation setup
setupObject = setupClass;
setupObject.saveObject.fileName = 'string';
setupObject.saveObject.saveData = false;
setupObject.totalTime = totalTime;
setupObject.totalTimeSteps = totalTime/timeStepSize;
setupObject.plotObject.flag = true;
setupObject.plotObject.steps = 0.2/timeStepSize;
setupObject.plotObject.keepFormerPlots = true;
setupObject.plotObject.plotInitialConfig = true;
setupObject.plotObject.postPlotType = 'disp';
setupObject.plotObject.stress.component = 1;
setupObject.plotObject.view = [0,90];
setupObject.plotObject.setXminXmax([-1;-2;0],[1;0;0]);
setupObject.newton.tolerance = newtonTolerance;
setupObject.integrator = integrator; 
dofObject = dofClass;   % required object for dof and object handling

%% Continuum object
stringObject = stringClass(dofObject,nDIM);
stringObject.elementDisplacementType = 'mixedPH'; %pH Formulation with right Cauchy-Green strains
stringObject.elementNameAdditionalSpecification = 'C';

% Material
stringObject.materialObject.name = material;
stringObject.materialObject.E = EA;
stringObject.materialObject.nu = nu;
stringObject.materialObject.lambda = EA*nu/((1+nu)*(1-2*nu));                 % Lame parameter
stringObject.materialObject.mu = EA/(2*(1+nu));
stringObject.materialObject.rho = rhoA; 
stringObject.numericalTangentObject.computeNumericalTangent = false;
stringObject.numericalTangentObject.showDifferences = false;

%% Spatial discretiaztion
endpoint = [sqrt(2)/2,-sqrt(2)/2];
length = norm(endpoint);
order  = 1;
number_of_gausspoints = 2;

% Mesh
[stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,endpoint);

% Shapefunctions
stringObject.dimension = 1;
stringObject.shapeFunctionObject.order = order;
stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;

%% Boundary conditions
% Dirichlet BC (fixed at x=0)
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = [1, 2];
dirichletObject.masterObject = stringObject; %% TO Do: use a term which is not racist

% neumann boundary conditions
nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = stringObject;
nodalLoadObject.loadVector = rhoA*[1; 1];
nodalLoadObject.timeFunction = str2func('@(t) sin(pi*t/0.2).*(t<=0.2+1e-12)'); 
nodalLoadObject.nodeList= size(stringObject.meshObject.nodes,1);

%Bodyforce
bodyForceObject = bodyForceClass(dofObject);
bodyForceObject.typeOfLoad = 'deadLoad';
bodyForceObject.masterObject = stringObject;
bodyForceObject.loadFunction = rhoA*[0; -9.81];
bodyForceObject.dimension = 2;
bodyForceObject.timeFunction = str2func('@(t) 1');
bodyForceObject.meshObject.edof = stringObject.meshObject.edof;
bodyForceObject.shapeFunctionObject.order = order;
bodyForceObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;

%% Solver
dofObject = runNewton(setupObject,dofObject);

%% Postprocessing

if postprocess

    % Plot the problem
    %plot(solidObject,setupObject)
    
    % Get energy quantities
    kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
    potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
    time = zeros(setupObject.totalTimeSteps+1,1);
    
    for j = 1:setupObject.totalTimeSteps+1
        time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
        kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
        potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).strainEnergy + bodyForceObject.ePot(j).externalEnergy;
    end
    
    % Compute increment
    figure()
    plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy)
    legend('T', 'V', 'H');
    matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'energy', '.tikz'], 'showInfo', false, 'floatformat', '%.7g')

    % Plot increment
    figure()
    energyIncrement = zeros(setupObject.totalTimeSteps,1);
    
    for i = 1:setupObject.totalTimeSteps
        energyIncrement(i) = abs((kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i))); 
    end
    plot(time(1:end-1),energyIncrement)
    legend('Energy increment');
    matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'Hdiff', '.tikz'], 'showInfo', false, 'floatformat', '%.4g')

end

