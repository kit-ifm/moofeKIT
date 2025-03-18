%% Visco-Elastic string, Left-to-Right-Maneuver
%
% problem: visco-hyperelastic or linear elastic string, maneuvering from
%          left to right using dirichlet and neumann BC, starts from static
%          equilibrium
% spatial discretization: linear Lagrange shapefunction for disp, 
%                         discont. constant shapefunction for strain,
%                         mixed formulation
% time discretization: endpoint, midpoint rule or discrete gradient
%
% user: Philipp Kinon

%% Main parameters
totalTime = 6;
newtonTolerance = 1e-8;
%integrator = 'Midpoint';
integrator = 'DiscreteGradient';
%material = 'StVenantVisco';
material = 'HyperelasticVisco';
EA = 20;
rhoA = 1;
tau_e = 0.02; % retardation time corresponding to elongation
etaA = tau_e*EA; % tau_e = eta_e/E
numberOfElementsOfCrime = 10;
load("equilibriumState.mat");

timeStepSize = 5e-2;
allDT = [5e-2, 5e-3, 5e-4, 2e-4];
u = cell(size(allDT));
error = zeros(size(allDT));
for k = 1:numel(allDT)

    timeStepSize = allDT(k);
    %% Simulation setup
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'string';
    setupObject.saveObject.saveData = false;
    setupObject.totalTime = totalTime;
    setupObject.totalTimeSteps = totalTime/timeStepSize;
    setupObject.plotObject.flag = false;
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
    stringObject.elementDisplacementType = 'mixedPH'; %pH Formulation with right Cauchy-Green strains
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
    endpoint = [0,-1];
    length_string = norm(endpoint);
    order  = 1;
    number_of_gausspoints = 2;
    
    % Mesh
    [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length_string,numberOfElementsOfCrime,order,startpoint,endpoint);
    stringObject.qN = nodesEquilibrium;
    dofObject.listContinuumObjects{1}.mixedFEObject.qR = mixedEquilibrium;

    % Shapefunctions
    stringObject.dimension = 1;
    stringObject.shapeFunctionObject.order = order;
    stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;
    
    %% Boundary conditions
    % Dirichlet BC (fixed at x=0)
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,2)==0);
    dirichletObject.nodalDof = [2];
    dirichletObject.masterObject = stringObject; %% TO Do: use a term which is not racist
    
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

    u{k} = stringObject.qN1;

end
      
%% Postprocessing
for l = 1:numel(allDT)
    error(l) = norm(u{l}-u{end})/norm(u{end});
end
figure()
loglog(allDT(1:end-1),error(1:end-1));
