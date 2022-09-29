% ROTATINGXELECTROTHERMOMECHANICS Script for preprocessing a dynamic
% coupled electro-thermo-mechanical simulation.
%
% REFERENCE
% In preparation
%
% SEE ALSO
% LShape
%
% CREATOR(S)
% Felix Zähringer, Marlon Franke

run('../../startUpMoofeKIT.m');

computeInitialConfiguration = true;

armWidth = 0.01;
armLength = 0.03;
bodyThickness = 0.002;

if computeInitialConfiguration
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.saveData = true;
    setupObject.saveObject.fileName = 'rotatingGripperBasis';
    lastTimeStep = 10;
    setupObject.saveObject.saveForTimeSteps = 1:lastTimeStep;
    setupObject.totalTimeSteps = lastTimeStep;
    setupObject.totalTime = 1;
    % setupObject.usePreconditioning = true;
    setupObject.newton.tolerance = 1e-5;
    setupObject.newton.maximumResiduumNorm = 1e30;
    setupObject.plotObject.flag = false;
    % setupObject.plotObject.postPlotType = 'phi';
    setupObject.plotObject.postPlotType = 'stress';
    % setupObject.plotObject.postPlotType = 'D1';
    % setupObject.plotObject.stress.component = -1;
    % setupObject.plotObject.view = [-36.632241204539113,25.552647559456744];
    % setupObject.integrator = 'Endpoint';
%     setupObject.integrator = 'Midpoint';
    setupObject.integrator = 'DiscreteGradient';

    dofObject = dofClass;   % required object for dof and object handling

    %% continuum Objects
    solidElectroThermoObject = solidElectroThermoClass(dofObject);
    [nodes, edof, bounEdofs] = meshGeneratorCube(2*armLength+armWidth, armWidth, bodyThickness, 14, 2, 2, 2, 1);
    nnodes=size(nodes(:,1),1);
    initialThermalField = 293.15*ones(nnodes,1);
    initialElectricalPotentialField = zeros(nnodes,1);
    initialNodes = [nodes, initialElectricalPotentialField, initialThermalField];
    solidElectroThermoObject.meshObject.nodes = initialNodes;
    solidElectroThermoObject.meshObject.edof = edof;

    solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
    solidElectroThermoObject.elementDisplacementType = 'mixedSC';
    solidElectroThermoObject.materialObject.rho = 1000;                                         % mass-density
    solidElectroThermoObject.materialObject.a = 25000;                                         % a
    solidElectroThermoObject.materialObject.b = 50000;                                         % b
    solidElectroThermoObject.materialObject.c1 = 500000;                                         % c
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
    solidElectroThermoObject.materialObject.thetaR = 273.15;                                  % Reference Temperature
    solidElectroThermoObject.materialObject.k0 = 0.23;                                          % thermal conductivity
    solidElectroThermoObject.materialObject.wk = 0;
    solidElectroThermoObject.materialObject.rhoSource = 0;
    solidElectroThermoObject.materialObject.RSource = 0;
    solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;

    solidElectroThermoObject.dimension = 3;
    solidElectroThermoObject.shapeFunctionObject.order = 2;
    solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidElectroThermoObject.mixedFEObject.condensation = true;
    solidElectroThermoObject.mixedFEObject.orderShapeFunction = 1;


    % initial velocity
    vN = zeros(size(solidElectroThermoObject.meshObject.nodes));
    omega = [0, 0, 2];
    for ii=1:size(vN, 1)
        vN(ii,1:3) = cross(solidElectroThermoObject.meshObject.nodes(ii,1:3), omega);
    end
    solidElectroThermoObject.vN = vN;

    % electrical dirichlet boundary conditions
    maxPotential = 5e4;
    tend = 10;
    tendLoad = 0.5;
    tendLoad2 = 0.5;
    tendLoad3 = 0.9;
    % timeFunction = @(t,Z) (maxPotential*sin(pi/2*t/tendLoad)).*(t >= 0).*(t <= tendLoad) + maxPotential.*(t > tendLoad).*(t <= tendLoad2) + (maxPotential*cos(0.5*pi*(t-tendLoad2)/(tendLoad3-tendLoad2))).*(t > tendLoad2).*(t <= tendLoad3);
    timeFunction = @(t,Z) (maxPotential*sin(pi/2*t/tendLoad)).*(t >= 0).*(t <= tendLoad) + maxPotential.*(t > tendLoad);
    nodeListDirichletObjects = cell(5, 1);

    % Mid
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 0);
    dirichletObject1.nodalDof = 4;
    dirichletObject1.masterObject = solidElectroThermoObject;
    nodeListDirichletObjects{1} = dirichletObject1.nodeList;

    % Top
    dirichletObject2 = dirichletClass(dofObject);
    dirichletObject2.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == (bodyThickness/2) & solidElectroThermoObject.meshObject.nodes(:,1) <= -(armWidth/2+1/3*armLength-1e-9));
    dirichletObject2.nodalDof = 4;
    dirichletObject2.timeFunction = timeFunction;
    dirichletObject2.masterObject = solidElectroThermoObject;
    nodeListDirichletObjects{2} = dirichletObject2.nodeList;

    dirichletObject3 = dirichletClass(dofObject);
    dirichletObject3.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == (bodyThickness/2) & solidElectroThermoObject.meshObject.nodes(:,1) >= (armWidth/2+1/3*armLength-1e-9));
    dirichletObject3.nodalDof = 4;
    dirichletObject3.timeFunction = timeFunction;
    dirichletObject3.masterObject = solidElectroThermoObject;
    nodeListDirichletObjects{3} = dirichletObject3.nodeList;

    % Bottom
    maxPotential = -maxPotential;
    timeFunction = @(t,Z) (maxPotential*sin(pi/2*t/tendLoad)).*(t >= 0).*(t <= tendLoad) + maxPotential.*(t > tendLoad);

    dirichletObject4 = dirichletClass(dofObject);
    dirichletObject4.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == -(bodyThickness/2) & solidElectroThermoObject.meshObject.nodes(:,1) <= -(armWidth/2-1e-9) & solidElectroThermoObject.meshObject.nodes(:,1) >= -(armWidth/2 + 2/3*armLength + 1e-9));
    dirichletObject4.nodalDof = 4;
    dirichletObject4.timeFunction = timeFunction;
    dirichletObject4.masterObject = solidElectroThermoObject;
    nodeListDirichletObjects{4} = dirichletObject4.nodeList;

    dirichletObject5 = dirichletClass(dofObject);
    dirichletObject5.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == -(bodyThickness/2) & solidElectroThermoObject.meshObject.nodes(:,1) >= (armWidth/2-1e-9) & solidElectroThermoObject.meshObject.nodes(:,1) <= (armWidth/2 + 2/3*armLength + 1e-9));
    dirichletObject5.nodalDof = 4;
    dirichletObject5.timeFunction = timeFunction;
    dirichletObject5.masterObject = solidElectroThermoObject;
    nodeListDirichletObjects{5} = dirichletObject5.nodeList;

    save('rotatingGripperInitialData.mat', 'nodeListDirichletObjects');

    %% solver
    warning off;
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    runNewton(setupObject,dofObject);
    
    clear;
end

%% Calculation convergence
computeConvergence = true;
if computeConvergence
    variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
    for jj=1:size(variableTimeStepSizeVector,1)
        variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
        % load Data of first steps
        load('rotatingGripperInitialData.mat');
        load('rotatingGripperBasis.mat', 'setupObject10', 'dofObject10');
        
        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = strcat('rotatingGripperDT', num2str(variableTimeStepSizeVector(jj)));
        setupObject.totalTime = 0.1;
        lastTimeStep = setupObject.totalTime/variableTimeStepSizeVector(jj);
        
        setupObject.totalTimeSteps = lastTimeStep;
        
        setupObject.saveObject.saveData = true;
        setupObject.saveObject.saveForTimeSteps = lastTimeStep;
        if mod(lastTimeStep, 1000) == 0
            setupObject.saveObject.saveForTimeSteps = 1000:1000:lastTimeStep;
        end
        setupObject.saveObject.exportLogFile = true;
        setupObject.newton.tolerance = 1e-5;
        setupObject.newton.maximumResiduumNorm = 1e30;
%         setupObject.usePreconditioning = true;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.makeMovie = false;
        % setupObject.plotObject.postPlotType = 'stress';
        % setupObject.plotObject.postPlotType = 'none';
        % setupObject.plotObject.stress.component = -1;
        setupObject.plotObject.view = [-36.632241204539113,25.552647559456744];
        % setupObject.integrator = 'Endpoint';
        % setupObject.integrator = 'Midpoint';
        setupObject.integrator = 'DiscreteGradient';
        
        if variableTimeStepSizeVector(jj) < 1e-4
            setupObject.integrator = 'Midpoint';
        else
            setupObject.integrator = 'DiscreteGradient';
        end
    
        dofObject = dofClass;   % required object for dof and object handling
    
        %% continuum Objects
        armWidth = 0.01;
        armLength = 0.03;
        bodyThickness = 0.002;

        solidElectroThermoObject = solidElectroThermoClass(dofObject);
        [nodes, edof, bounEdofs] = meshGeneratorCube(2*armLength+armWidth, armWidth, bodyThickness, 14, 2, 2, 2, 1);
        nnodes=size(nodes(:,1),1);
        initialThermalField = 293.15*ones(nnodes,1);
        initialElectricalPotentialField = zeros(nnodes,1);
        initialNodes = [nodes, initialElectricalPotentialField, initialThermalField];
        solidElectroThermoObject.meshObject.nodes = initialNodes;
        solidElectroThermoObject.meshObject.edof = edof;

        solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
        solidElectroThermoObject.elementDisplacementType = 'mixedSC';
        solidElectroThermoObject.materialObject.rho = 1000;                                         % mass-density
        solidElectroThermoObject.materialObject.a = 25000;                                         % a
        solidElectroThermoObject.materialObject.b = 50000;                                         % b
        solidElectroThermoObject.materialObject.c1 = 500000;                                         % c
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
        solidElectroThermoObject.materialObject.thetaR = 273.15;                                  % Reference Temperature
        solidElectroThermoObject.materialObject.k0 = 0.23;                                          % thermal conductivity
        solidElectroThermoObject.materialObject.wk = 0;
        solidElectroThermoObject.materialObject.rhoSource = 0;
        solidElectroThermoObject.materialObject.RSource = 0;
        solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
        
        solidElectroThermoObject.dimension = 3;
        solidElectroThermoObject.shapeFunctionObject.order = 2;
        solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
        solidElectroThermoObject.mixedFEObject.condensation = true;
        solidElectroThermoObject.mixedFEObject.orderShapeFunction = 1;
    
        solidElectroThermoObject.mixedFEObject.qN = dofObject10.listContinuumObjects{1}.mixedFEObject.qN1;
        solidElectroThermoObject.qN = dofObject10.listContinuumObjects{1}.qN1;
        solidElectroThermoObject.vN = dofObject10.listContinuumObjects{1}.vN1;

        % electrical dirichlet boundary conditions
        maxPotential = 5e4;
        timeFunction = @(t,Z) maxPotential;
    
        % Mid
        dirichletObject1 = dirichletClass(dofObject);
        dirichletObject1.nodeList = nodeListDirichletObjects{1};
        dirichletObject1.nodalDof = 4;
        dirichletObject1.masterObject = solidElectroThermoObject;
    
        % Top
        dirichletObject2 = dirichletClass(dofObject);
        dirichletObject2.nodeList = nodeListDirichletObjects{2};
        dirichletObject2.nodalDof = 4;
        dirichletObject2.timeFunction = timeFunction;
        dirichletObject2.masterObject = solidElectroThermoObject;
    
        dirichletObject3 = dirichletClass(dofObject);
        dirichletObject3.nodeList = nodeListDirichletObjects{3};
        dirichletObject3.nodalDof = 4;
        dirichletObject3.timeFunction = timeFunction;
        dirichletObject3.masterObject = solidElectroThermoObject;
    
        % Bottom
        maxPotential = -maxPotential;
        timeFunction = @(t,Z) maxPotential;
    
        dirichletObject4 = dirichletClass(dofObject);
        dirichletObject4.nodeList = nodeListDirichletObjects{4};
        dirichletObject4.nodalDof = 4;
        dirichletObject4.timeFunction = timeFunction;
        dirichletObject4.masterObject = solidElectroThermoObject;
    
        dirichletObject5 = dirichletClass(dofObject);
        dirichletObject5.nodeList = nodeListDirichletObjects{5};
        dirichletObject5.nodalDof = 4;
        dirichletObject5.timeFunction = timeFunction;
        dirichletObject5.masterObject = solidElectroThermoObject;

        %% solver
        warning off;
        parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
        parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
        runNewton(setupObject,dofObject);
        
        clear;
    end
end