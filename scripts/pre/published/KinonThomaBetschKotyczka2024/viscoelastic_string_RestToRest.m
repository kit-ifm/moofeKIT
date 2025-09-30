%% Visco-Elastic string, Left-to-Right-Maneuver
%
% problem: visco-hyperelastic string, maneuvering from
%          left to right using dirichlet and neumann BC, starts from static
%          equilibrium
% spatial discretization: linear Lagrange shapefunction for disp,
%                         discont. constant shapefunction for strain,
%                         mixed formulation
% time discretization: midpoint rule or discrete gradient
%
% user: Philipp Kinon

%% Main parameters
totalTime = 6;
timeStepSize = 5e-2;
newtonTolerance = 1e-10;
%integrator = 'Midpoint';
integrator = 'DiscreteGradient';
material = 'Hyperelastic';
EA = 20;
rhoA = 1;
tau_e = 0.02; % retardation time corresponding to elongation
etaA = tau_e*EA; % tau_e = eta_e/E
numberOfElementsOfCrime = 10;
postprocess = true;
new_calculation = true;
exportFolder = '../../../../generatedFiles/imageFiles/';
load("equilibriumState.mat");

if new_calculation

    %% Simulation setup
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'string';
    setupObject.saveObject.saveData = false;
    setupObject.totalTime = totalTime;
    setupObject.totalTimeSteps = totalTime/timeStepSize;
    setupObject.plotObject.flag = true;
    setupObject.plotObject.steps = 0.4/timeStepSize;
    setupObject.plotObject.keepFormerPlots = true;
    setupObject.plotObject.plotInitialConfig = true;
    setupObject.plotObject.postPlotType = 'disp';
    setupObject.plotObject.stress.component = 1;
    setupObject.plotObject.view = [0,90];
    setupObject.plotObject.setXminXmax([-0.1;-2;0],[2;0;0]);
    setupObject.newton.tolerance = newtonTolerance;
    setupObject.integrator = integrator;
    dofObject = dofClass;
    %% Continuum object
    stringObject = stringClass(dofObject,2);
    stringObject.elementDisplacementType = 'mixedPHViscoPT'; %pH Formulation with right Cauchy-Green strains
    stringObject.elementNameAdditionalSpecification = 'C';

    % Material
    stringObject.materialObject.name = material;
    stringObject.materialObject.EA = EA;
    stringObject.materialObject.etaA = etaA;
    stringObject.materialObject.rho = rhoA;
    stringObject.numericalTangentObject.computeNumericalTangent = false;
    stringObject.numericalTangentObject.showDifferences = false;

    %% Spatial discretiaztion
    endpoint = [0,-1];
    length = norm(endpoint);
    order  = 1;
    number_of_gausspoints = 2;

    % Mesh
    [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,[0,0],endpoint);
    stringObject.qN = nodesEquilibrium;
    dofObject.listContinuumObjects{1}.mixedFEObject.qN = mixedEquilibrium;

    % Shapefunctions
    stringObject.dimension = 1;
    stringObject.shapeFunctionObject.order = order;
    stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;

    %% Boundary conditions
    % Dirichlet BC (fixed at x=0)
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,2)==0);
    dirichletObject.nodalDof = [2];
    dirichletObject.masterObject = stringObject;

    % neumann boundary conditions
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = stringObject;
    nodalLoadObject.loadVector = rhoA*[1; 0];
    nodalLoadObject.timeFunction = str2func('@(t) sin(pi*(t)/2).*(t<=4+1e-12)');
    nodalLoadObject.nodeList= 1;

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

end

%% Postprocessing

if postprocess

    % Plot the problem
    plot(stringObject,setupObject)
    %matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'snapshots', '.tikz'], 'showInfo', false, 'floatformat', '%.4g')


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
        potentialEnergy(j) = dofObject.listContinuumObjects{1,1}.elementData(j).strainEnergy + bodyForceObject.elementData(j).externalEnergy;
        dissipatedEnergy(j) = dofObject.listContinuumObjects{1,1}.elementData(j).dissipatedEnergy;
    end
    for i = 1:setupObject.totalTimeSteps
        energyIncrement(i) = (kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i));
        totalDissipatedEnergy(i+1) = totalDissipatedEnergy(i) + dissipatedEnergy(i+1);
    end

    figure()
    plot(time(1:end-1),dissipatedEnergy(2:end))
    legend('W_{diss}^n');
    matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'Wdiss', '.tikz'], 'showInfo', false, 'floatformat', '%.4g')


    figure()
    plot(time, kineticEnergy+potentialEnergy, time, totalDissipatedEnergy, time, kineticEnergy+potentialEnergy+totalDissipatedEnergy)
    legend('H', 'E_{diss}','H+E_{diss}');
    matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'H_Ediss', '.tikz'], 'showInfo', false, 'floatformat', '%.4g')


    figure()
    semilogy(time(1:end-1),abs(energyIncrement + dissipatedEnergy(2:end)))
    legend('H^{n+1}-H^n + W_{diss}^n');
    grid on
    matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [exportFolder, 'H_W_diff', '.tikz'], 'showInfo', false, 'floatformat', '%.4g')


end