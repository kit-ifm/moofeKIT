%% Visco-Elastic string pendulum
%
% problem: visco-hyperelastic or linear elastic string, fixed on one end, 
%          time-dependent force on the other
% spatial discretization: linear Lagrange shapefunction for disp
%                         discont. constant ansatz for strain, 
%                         mixed formulation
% time discretization: endpoint, midpoint rule or discrete gradient
%
% user: Philipp Kinon

%% Main parameters
totalTime = 1;
timeStepSize = 1e-2;
newtonTolerance = 1e-9;
%integrator = 'Midpoint';
integrator = 'DiscreteGradient';
%material = 'StVenantVisco';
%material = 'HyperelasticVisco';
material = 'Hyperelastic';
EA = 20;
rhoA = 1;
tau_e = 0.02; % retardation time corresponding to elongation
etaA = tau_e*EA; % tau_e = eta_e/E
numberOfElementsOfCrime = 10;
postprocess = true;

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
dofObject = dofClass;
%% Continuum object
stringObject = stringClass(dofObject,2);
%stringObject.elementDisplacementType = 'mixedPH'; %pH Formulation with right Cauchy-Green strains
stringObject.elementDisplacementType = 'mixedPHViscoPT'; %pH Formulation with right Cauchy-Green strains and Poynting-Thompson viscoelasticity
stringObject.elementNameAdditionalSpecification = 'C';

% Material
stringObject.materialObject.name = material;
stringObject.materialObject.EA = EA;
stringObject.materialObject.etaA = etaA;
stringObject.materialObject.rho = rhoA; 
stringObject.numericalTangentObject.computeNumericalTangent = false;
stringObject.numericalTangentObject.showDifferences = false;

%% Spatial discretiaztion
startpoint = [0,0];
endpoint = [sqrt(2)/2,-sqrt(2)/2];
length = norm(endpoint);
order  = 1;
number_of_gausspoints = 2;

% Mesh
[stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,startpoint,endpoint);

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
    plot(stringObject,setupObject)
    
    % Get energy quantities
    kineticEnergy = zeros(setupObject.totalTimeSteps+1,1);
    potentialEnergy = zeros(setupObject.totalTimeSteps+1,1);
    dissipatedEnergy = zeros(setupObject.totalTimeSteps+1,1);
    time = zeros(setupObject.totalTimeSteps+1,1);
    totalDissipatedEnergy = zeros(setupObject.totalTimeSteps+1,1);
    energyIncrement = zeros(setupObject.totalTimeSteps,1);
       
    for j = 1:setupObject.totalTimeSteps+1
        time(j) = (j-1)*setupObject.totalTime/setupObject.totalTimeSteps;
        kineticEnergy(j) = dofObject.postDataObject.energyJournal(j).EKin;
        potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).strainEnergy + bodyForceObject.ePot(j).externalEnergy;
        dissipatedEnergy(j) = dofObject.listContinuumObjects{1,1}.ePot(j).dissipatedEnergy;
    end
    for i = 1:setupObject.totalTimeSteps
        energyIncrement(i) = (kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i)); 
        totalDissipatedEnergy(i+1) = totalDissipatedEnergy(i) + dissipatedEnergy(i+1);
    end
    
    
    % Compute increment
    figure()
    plot(time, kineticEnergy, time, potentialEnergy, time, kineticEnergy+potentialEnergy, time, totalDissipatedEnergy, time, kineticEnergy+potentialEnergy+totalDissipatedEnergy)
    legend('T', 'V', 'H', 'E_{diss}','H+E_{diss}');
    
    % Plot increment
    figure()
    plot(time(1:end-1),energyIncrement,time(1:end-1),dissipatedEnergy(2:end))
    legend('H^{n+1}-H^n','W_{diss}^n');

    figure()
    plot(time(1:end-1),energyIncrement + dissipatedEnergy(2:end))
    legend('H^{n+1}-H^n + W_{diss}^n');


end