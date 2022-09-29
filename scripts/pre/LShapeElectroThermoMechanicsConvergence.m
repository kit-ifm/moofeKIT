% LSHAPEELECTROTHERMOMECHANICS Script for preprocessing a dynamic coupled 
% electro-thermo-mechanical simulation.
% 
% FORMULATION
% Mixed finite element formulation for nonlinear electro-thermo-mechanical 
% processes is pursued. This is accomplished under the assumption of a 
% hyperelastic, isotropic material and based on a polyconvex inspired 
% internal energy function. The discretization in time is pursued with an 
% energy and momentum conserving scheme. 
% 
% 
% REFERENCE
% https://doi.org/10.1016/j.cma.2021.114298
% https://doi.org/10.1007/BF00913408
% 
% SEE ALSO
% LShape, 
% cooksMembrane
% 
% CREATOR(S) 
% Marlon Franke

run('../../startUpMoofeKIT.m');

computeInitialConfiguration = true;

if computeInitialConfiguration
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeElectroThermoConvergenceBasis';
setupObject.saveObject.saveData = true;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 5;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = false;
% setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.postPlotType = 'phi';
setupObject.plotObject.postPlotType = 'stress';
% setupObject.plotObject.postPlotType = 'D1';
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityCoarse.inp');
solidElectroThermoObject = solidElectroThermoClass(dofObject);
solidElectroThermoObject.meshObject.nodes = abaqusMeshData.QNODE;
nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);
nodelistoben = abaqusMeshData.SUBSETS(3).EDOF(:);
nodelistunten = abaqusMeshData.SUBSETS(4).EDOF(:);
sizenodesoben = size(nodelistoben(:,1),1);
for l = 1:sizenodesoben
    initialThermalField(nodelistoben(l))=300;
    initialThermalField(nodelistunten(l))=250;
end

initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.meshObject.edof = abaqusMeshData.EDOF;
solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.rho = 100;                                          % mass-density
solidElectroThermoObject.materialObject.a = 831.25;                                         % a
solidElectroThermoObject.materialObject.b = 166.25;                                         % b
solidElectroThermoObject.materialObject.c1 = 10000;                                         % c
solidElectroThermoObject.materialObject.c2 = 5209;                                          % e
solidElectroThermoObject.materialObject.d1 = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12);                             	% epsilon_0
solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 1500;                                        % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                             % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 293.15;                                  % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 0.23;                                          % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.orderShapeFunction = 1;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = false;

bcTimeEnd = 5;
FA = [256;512;768];

neumannObject1 = neumannClass(dofObject);
neumannObject1.field = 'mechanical';
neumannObject1.typeOfLoad = 'deadLoad';
neumannObject1.masterObject = solidElectroThermoObject;
neumannObject1.forceVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.SUBSETS(3).EDOF;

neumannObject2 = neumannClass(dofObject);
neumannObject2.field = 'mechanical';
neumannObject2.typeOfLoad = 'deadLoad';
neumannObject2.masterObject = solidElectroThermoObject;
neumannObject2.forceVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.SUBSETS(4).EDOF;

% electrical dirichlet boundary conditions
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 0);
nodesDirichlet1 = dirichletObject1.nodeList;
dirichletObject1.nodalDof = 4;
dirichletObject1.masterObject = solidElectroThermoObject;
dirichletObject1.timeFunction = str2func('@(t,Z) (6e6).*(t > 5) + (6e6*sin(pi/2*t/5)).*(t >= 0).*(t <= 5)');  

dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 1.5);
nodesDirichlet2 = dirichletObject2.nodeList;
dirichletObject2.nodalDof = 4;
dirichletObject2.masterObject = solidElectroThermoObject;

save('LShapeElectroThermoConvergenceNodeLists', 'nodesDirichlet1', 'nodesDirichlet2');

%% solver
warning off;
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix');
dofObject = runNewton(setupObject,dofObject);

end

clear;

variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
for jj=1:size(variableTimeStepSizeVector,1)
    variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
    % load Data of first steps
    load('LShapeElectroThermoConvergenceBasis', 'setupObject10', 'dofObject10');
    
    setupObject = setupClass;
    setupObject.saveObject.fileName = strcat('LShapeElectroThermoConvergence', num2str(variableTimeStepSizeVector(jj)));
    setupObject.totalTime = 0.1;
    lastTimeStep = setupObject.totalTime/variableTimeStepSizeVector(jj);
    setupObject.totalTimeSteps = lastTimeStep;
    
    setupObject.saveObject.saveData = true;
    setupObject.saveObject.saveForTimeSteps = lastTimeStep;
    if mod(lastTimeStep, 1000) == 0
        setupObject.saveObject.saveForTimeSteps = 1000:1000:lastTimeStep;
    end
    setupObject.plotObject.flag = false;
    % setupObject.plotObject.stress.component = -1;
    setupObject.plotObject.view = [-0.4,90];
    setupObject.newton.tolerance = 1.2e-3;
    % setupObject.integrator = 'Endpoint';
    if variableTimeStepSizeVector(jj) <= 1e-4
        setupObject.integrator = 'Midpoint';
    else
        setupObject.integrator = 'DiscreteGradient';
    end

    dofObject = dofClass;   % required object for dof and object handling
    
    %% continuum Objects
    % abaqus mesh
    abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityCoarse.inp');
    solidElectroThermoObject = solidElectroThermoClass(dofObject);
    solidElectroThermoObject.meshObject.nodes = abaqusMeshData.QNODE;
    nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
    initialThermalField = 293.15*ones(nnodes,1);
    nodelistoben = abaqusMeshData.SUBSETS(3).EDOF(:);
    nodelistunten = abaqusMeshData.SUBSETS(4).EDOF(:);
    sizenodesoben = size(nodelistoben(:,1),1);
    for l = 1:sizenodesoben
        initialThermalField(nodelistoben(l))=300;
        initialThermalField(nodelistunten(l))=250;
    end

    initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
    solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
    solidElectroThermoObject.meshObject.edof = abaqusMeshData.EDOF;
    solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
    solidElectroThermoObject.elementDisplacementType = 'mixedSC';
    solidElectroThermoObject.materialObject.rho = 100;                                          % mass-density
    solidElectroThermoObject.materialObject.a = 831.25;                                         % a
    solidElectroThermoObject.materialObject.b = 166.25;                                         % b
    solidElectroThermoObject.materialObject.c1 = 10000;                                         % c
    solidElectroThermoObject.materialObject.c2 = 5209;                                          % e
    solidElectroThermoObject.materialObject.d1 = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
    solidElectroThermoObject.materialObject.d2 = 0;
    solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12);                             	% epsilon_0
    solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
    solidElectroThermoObject.materialObject.e2 = 0;
    solidElectroThermoObject.materialObject.ee = 0;
    solidElectroThermoObject.materialObject.kappa = 1500;                                        % heat capacity
    solidElectroThermoObject.materialObject.strongcomp = false;
    solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                             % coupling parameter
    solidElectroThermoObject.materialObject.thetaR = 293.15;                                  % Reference Temperature
    solidElectroThermoObject.materialObject.k0 = 0.23;                                          % thermal conductivity
    solidElectroThermoObject.materialObject.wk = 0;
    solidElectroThermoObject.materialObject.rhoSource = 0;
    solidElectroThermoObject.materialObject.RSource = 0;
    solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
    solidElectroThermoObject.dimension = 3;
    solidElectroThermoObject.shapeFunctionObject.order = 2;
    solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidElectroThermoObject.mixedFEObject.condensation = true;
    solidElectroThermoObject.mixedFEObject.orderShapeFunction = 1;
    solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
    solidElectroThermoObject.numericalTangentObject.showDifferences = false;
    
    solidElectroThermoObject.vN = dofObject10.listContinuumObjects{1}.vN1;
    solidElectroThermoObject.qN = dofObject10.listContinuumObjects{1}.qN1;
    solidElectroThermoObject.mixedFEObject.qN = dofObject10.listContinuumObjects{1}.mixedFEObject.qN1;


    % electrical dirichlet boundary conditions
    load('LShapeElectroThermoConvergenceNodeLists');
    
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.nodeList = nodesDirichlet1;
    dirichletObject1.nodalDof = 4;
    dirichletObject1.masterObject = solidElectroThermoObject;
    dirichletObject1.timeFunction = str2func('@(t,Z) (6e6)');  

    dirichletObject2 = dirichletClass(dofObject);
    dirichletObject2.nodeList = nodesDirichlet2;
    dirichletObject2.nodalDof = 4;
    dirichletObject2.masterObject = solidElectroThermoObject;
    
    

    %% solver
    warning off;
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix');
    dofObject = runNewton(setupObject,dofObject);

    clear;
end
