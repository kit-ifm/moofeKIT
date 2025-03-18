disp('=======================================');
disp('Now testing: Beam scripts');
disp('=======================================');

% turn off figures
set(0,'DefaultFigureVisible','off');

% name test cases
testCases = {'static_bernoulli_standard','static_timoshenko_standard','static_timoshenko_PH',...
    'static_timoshenko_standard_purebending','static_timoshenko_PH_purebending', ...
    'static_timoshenko_standard_forceAndMoment','static_timoshenko_PH_forceAndMoment', ...
    'dynamic_bernoulli_standard','dynamic_timoshenko_standard','dynamic_timoshenko_PH', ...
    'static_timoshenko_ANS', 'static_timoshenko_PetrovGalerkinANS'};

for i = 1:numel(testCases)
    % Catch data for test case
    materialName = 'Hooke';
    elementNameAdditionalSpecification = '';
    switch testCases{i}
        case 'static_bernoulli_standard'
            widthToLength = 1/50;
            beamType = 'Bernoulli';
            formulationType = 'displacement';
            nGP = 2;
            timeSteps = 1;
            Q = [100;0];
            q = [1;0];
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 16.5344;
            error_tolerance = 1e-4;
        case 'static_timoshenko_standard'
            widthToLength = 1/20;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            nGP = 1;
            timeSteps = 1;
            Q = [100;0];
            q = [1;0];
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 0.411324867724855;
            error_tolerance = 1e-12;
        case 'static_timoshenko_PH'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'mixedPH';
            nGP = 2;
            timeSteps = 1;
            Q = -[100;0]; % positive z-direction is upwards
            q = [-1;0];
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = -2.2457e-04;
            error_tolerance = 1e-8;
        case 'static_timoshenko_standard_purebending'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            nGP = 1;
            timeSteps = 1;
            Q = [0;0];
            q = [0;-10]; % positive phi-direction is counterclockwise
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 1.5e-05; % positive z-direction is downwards
            error_tolerance = 1e-9;
        case 'static_timoshenko_PH_purebending'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'mixedPH';
            nGP = 2;
            timeSteps = 1;
            Q = [0;0];
            q = [0;-10]; % positive phi-direction is counterclockwise
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = -1.5e-05; % positive z-direction is upwards
            error_tolerance = 1e-9;
        case 'static_timoshenko_standard_forceAndMoment'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            nGP = 1;
            timeSteps = 1;
            Q = [100;0];
            q = [1;-10]; % positive phi-direction is counterclockwise
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 2.3957e-04; % positive z-direction is downwards
            error_tolerance = 1e-8;
        case 'static_timoshenko_PH_forceAndMoment'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'mixedPH';
            nGP = 2;
            timeSteps = 1;
            Q = -[100;0];
            q = [-1;-10]; % positive phi-direction is counterclockwise
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = -2.3957e-04; % positive z-direction is upwards
            error_tolerance = 1e-8;
        case 'dynamic_bernoulli_standard'
            widthToLength = 1/50;
            beamType = 'Bernoulli';
            formulationType = 'displacement';
            nGP = 2;
            timeSteps = 10;
            Q = [1000;0];
            q = [10;20];
            rhoA = 1;
            integrator = 'Midpoint';
            timeFunction = str2func('@(t) sin(pi*t/0.5).*(t<=0.5+1e-12)');
            correct_value = [];
            error_tolerance = [];
        case 'dynamic_timoshenko_standard'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            nGP = 1;
            timeSteps = 10;
            Q = [1000;0];
            q = [10;20];
            rhoA = 1;
            integrator = 'Midpoint';
            timeFunction = str2func('@(t) sin(pi*t/0.5).*(t<=0.5+1e-12)');
            correct_value = [];
            error_tolerance = [];
        case 'dynamic_timoshenko_PH'
            widthToLength = 1/3;
            beamType = 'Timoshenko';
            formulationType = 'mixedPH';
            nGP = 2;
            timeSteps = 10;
            Q = [1000;0];
            q = [10;20];
            rhoA = 1;
            integrator = 'Midpoint';
            timeFunction = str2func('@(t) sin(pi*t/0.5).*(t<=0.5+1e-12)');
            correct_value = [];
            error_tolerance = [];
        case 'static_timoshenko_ANS'
            widthToLength = 1/20;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            elementNameAdditionalSpecification = 'ANS';
            nGP = 2;
            timeSteps = 1;
            Q = [100;0];
            q = [1;0];
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 0.411324867724855;
            error_tolerance = 1e-12;
        case 'static_timoshenko_PetrovGalerkinANS'
            widthToLength = 1/20;
            beamType = 'Timoshenko';
            formulationType = 'displacement';
            elementNameAdditionalSpecification = 'PetrovGalerkinANS';
            nGP = 2;
            timeSteps = 1;
            Q = [100;0];
            q = [1;0];
            rhoA = 0;
            integrator = 'Endpoint';
            timeFunction = str2func('@(t) 1');
            correct_value = 0.411324867724855;
            error_tolerance = 1e-12;
    end

    %% setup (setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'testbeams';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = timeSteps;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.newton.tolerance = 1e-9;
    setupObject.integrator = integrator;
    setupObject.plotObject.flag = false;

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    beamObject = beamClass(dofObject);
    beamObject.materialObject.name = materialName;
    beamObject.theory = beamType;
    beamObject.elementDisplacementType = formulationType;
    beamObject.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
    lengthBeam = 100;
    widthBeam = widthToLength * lengthBeam;
    beamObject.materialObject.E = 2.1 * 10^6;
    nu = 0.3;
    beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
    beamObject.materialObject.I = widthBeam^4 / 12;
    beamObject.materialObject.A = widthBeam^2;
    beamObject.materialObject.rho = rhoA;
    
    numberOfElements = 3;
    [nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
    nodes = nodes + 1 / 2 * lengthBeam;
    beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];
    
    beamObject.dimension = 1;
    beamObject.shapeFunctionObject.order = 1;
    beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
    
    % dirichlet boundaries
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = [1, 2];
    dirichletObject.masterObject = beamObject;
    
    % neumann boundary
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = beamObject;
    nodalLoadObject.loadVector = Q;
    nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);
    nodalLoadObject.timeFunction = timeFunction;

    % distributed load
    distributedForceObject = bodyForceClass(dofObject);
    distributedForceObject.typeOfLoad = 'deadLoad';
    distributedForceObject.masterObject = beamObject;
    distributedForceObject.loadFunction = q;
    distributedForceObject.timeFunction = timeFunction;
    distributedForceObject.dimension = 2; % second dimension required to fix bodyforce vector
    distributedForceObject.meshObject.edof = beamObject.meshObject.edof;
    distributedForceObject.shapeFunctionObject.order = 1;
    distributedForceObject.shapeFunctionObject.numberOfGausspoints = nGP;

    %% solver
    dofObject = runNewton(setupObject, dofObject);
    
    %% Check result
    if ~isempty(correct_value)
        displacementLastNode = dofObject.listContinuumObjects{1}.qN1(end, 1);
        assert(abs(displacementLastNode - correct_value) <= error_tolerance, ['Wrong results for testcase: ',testCases{i}]);
    end
    
end

close all