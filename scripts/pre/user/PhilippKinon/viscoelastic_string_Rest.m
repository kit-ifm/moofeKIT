%% Visco-Elastic string, static
%
% problem: visco-hyperelastic or linear elastic string, fixed on top,
%          static case
% spatial discretization: linear Lagrange shapefunction for disp
%                         discont. constant ansatz for strain, 
%                         mixed formulation
% time discretization: endpoint, midpoint rule or discrete gradient
%
% user: Philipp Kinon

%% Main parameters
totalTime = 10;
timeStepSize = 1e-1;
newtonTolerance = 1e-10;
%integrator = 'Midpoint';
%integrator = 'DiscreteGradient';
integrator = 'Endpoint';
%material = 'StVenantVisco';
material = 'HyperelasticVisco';
EA = 20;
rhoA = 1;
tau_e = 0.02; % retardation time corresponding to elongation
etaA = tau_e*EA; % tau_e = eta_e/E
numberOfElementsOfCrime = 10;
new_calculation = true;
exportFolder = '../../../../generatedFiles/imageFiles/';

if new_calculation

    %% Simulation setup
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'string';
    setupObject.saveObject.saveData = false;
    setupObject.totalTime = totalTime;
    setupObject.totalTimeSteps = totalTime/timeStepSize;
    setupObject.plotObject.flag = true;
    setupObject.plotObject.steps = setupObject.totalTimeSteps;
    setupObject.plotObject.keepFormerPlots = false;
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
    stringObject.elementDisplacementType = 'mixedPH'; %pH Formulation with right Cauchy-Green strains
    stringObject.elementNameAdditionalSpecification = 'C';
    
    % Material
    stringObject.materialObject.name = material;
    stringObject.materialObject.EA = EA;
    stringObject.materialObject.etaA = etaA;
    stringObject.materialObject.rho = 0; 
    stringObject.numericalTangentObject.computeNumericalTangent = false;
    stringObject.numericalTangentObject.showDifferences = false;
    
    %% Spatial discretiaztion
    startpoint = [0,0];
    endpoint = [0,-1];
    length = norm(endpoint);
    order  = 1;
    number_of_gausspoints = 2;
    
    % Mesh
    [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,startpoint,endpoint);
    
    % Specify mixed quantities different from 1 (required for static case)
    dofObject.listContinuumObjects{1}.mixedFEObject.qR = 1.5*ones(numberOfElementsOfCrime,1);

    % Shapefunctions
    stringObject.dimension = 1;
    stringObject.shapeFunctionObject.order = order;
    stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;
    
    %% Boundary conditions
    % Dirichlet BC (fixed at x=0)
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,2)==0);
    dirichletObject.nodalDof = [1 2];
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

end

nodesEquilibrium = reshape(dofObject.qN1(1:2*(numberOfElementsOfCrime+1)),[2,numberOfElementsOfCrime+1])';
mixedEquilibrium = dofObject.qN1(2*(numberOfElementsOfCrime+1)+1:end);
save("equilibriumState.mat","nodesEquilibrium","mixedEquilibrium")