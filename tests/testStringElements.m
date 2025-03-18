% unit test for string elements

disp('=======================================');
disp('Now testing: string elements');
disp('=======================================');

testCases = {'elastic_dynamic','visco_static','visco_dynamics'};

for ii = 1:size(testCases, 2)
    disp('---------------------------------------');
    disp(['- Current Testcase: ', testCases{ii}]);
    disp('---------------------------------------');

    if strcmp(testCases{ii},'elastic_dynamic')
        %% Main parameters
        totalTime = 1;
        timeStepSize = 1e-2;
        newtonTolerance = 1e-7;
        integrator = 'DiscreteGradient';
        material = 'Hyperelastic';
        EA = 10;
        rhoA = 1;
        numberOfElementsOfCrime = 30;
        nDIM = 2;
        
        %% Simulation setup
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'string';
        setupObject.saveObject.saveData = false;
        setupObject.totalTime = totalTime;
        setupObject.totalTimeSteps = totalTime/timeStepSize;
        setupObject.plotObject.flag = false;
        setupObject.newton.tolerance = newtonTolerance;
        setupObject.integrator = integrator; 
        dofObject = dofClass;   
        
        %% Continuum object
        stringObject = stringClass(dofObject,nDIM);
        stringObject.elementDisplacementType = 'mixedPH'; 
        stringObject.elementNameAdditionalSpecification = 'C';
        
        % Material
        stringObject.materialObject.name = material;
        stringObject.materialObject.EA = EA;
        stringObject.materialObject.rho = rhoA; 
        stringObject.numericalTangentObject.computeNumericalTangent = false;
        stringObject.numericalTangentObject.showDifferences = false;
        
        %% Spatial discretiaztion
        endpoint = [sqrt(2)/2,-sqrt(2)/2];
        length = norm(endpoint);
        order  = 1;
        number_of_gausspoints = 2;
        
        % Mesh
        [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,[0;0],endpoint);
        
        % Shapefunctions
        stringObject.dimension = 1;
        stringObject.shapeFunctionObject.order = order;
        stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;
        
        %% Boundary conditions
        dirichletObject = dirichletClass(dofObject);
        dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,1)==0);
        dirichletObject.nodalDof = [1, 2];
        dirichletObject.masterObject = stringObject; 
        
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
    
    elseif strcmp(testCases{ii},'visco_static')
    
        %% Main parameters
        totalTime = 1;
        timeStepSize = 1e-1;
        newtonTolerance = 1e-10;
        integrator = 'Endpoint';
        material = 'HyperelasticVisco';
        EA = 20;
        rhoA = 1;
        tau_e = 0.02; % retardation time corresponding to elongation
        etaA = 0.02*EA; % tau_e = eta_e/E
        numberOfElementsOfCrime = 10;
    
        %% Simulation setup
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'string';
        setupObject.saveObject.saveData = false;
        setupObject.totalTime = totalTime;
        setupObject.totalTimeSteps = totalTime/timeStepSize;
        setupObject.plotObject.flag = false;
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
        endpoint = [0,-1];
        length = norm(endpoint);
        order  = 1;
        number_of_gausspoints = 2;
        
        % Mesh
        [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,[0;0],endpoint);
        
        % Specify mixed quantities different from 1 (required for static case)
        dofObject.listContinuumObjects{1}.mixedFEObject.qR = 1.5*ones(numberOfElementsOfCrime,1);
    
        % Shapefunctions
        stringObject.dimension = 1;
        stringObject.shapeFunctionObject.order = order;
        stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;
        
        %% Boundary conditions
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
    
    elseif strcmp(testCases{ii},'visco_dynamics')
    
        %% Main parameters
        totalTime = 1;
        timeStepSize = 1e-2;
        newtonTolerance = 1e-9;
        integrator = 'DiscreteGradient';
        material = 'HyperelasticVisco';
        EA = 20;
        rhoA = 1;
        tau_e = 0.02; % retardation time corresponding to elongation
        etaA = tau_e*EA; % tau_e = eta_e/E
        numberOfElementsOfCrime = 10;
        
        %% Simulation setup
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'string';
        setupObject.saveObject.saveData = false;
        setupObject.totalTime = totalTime;
        setupObject.totalTimeSteps = totalTime/timeStepSize;
        setupObject.plotObject.flag = false;
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
        endpoint = [sqrt(2)/2,-sqrt(2)/2];
        length = norm(endpoint);
        order  = 1;
        number_of_gausspoints = 2;
        
        % Mesh
        [stringObject.meshObject.nodes,stringObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElementsOfCrime,order,[0;0],endpoint);
        
        % Shapefunctions
        stringObject.dimension = 1;
        stringObject.shapeFunctionObject.order = order;
        stringObject.shapeFunctionObject.numberOfGausspoints = number_of_gausspoints;
        
        %% Boundary conditions
        dirichletObject = dirichletClass(dofObject);
        dirichletObject.nodeList = find(stringObject.meshObject.nodes(:,1)==0);
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

end