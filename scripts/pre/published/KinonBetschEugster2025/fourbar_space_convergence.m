%% Dynamical GE beam multibody system
clearvars; close all;

%% Parameters
order  = 1;
beamType = 'GeometricallyExact';
nGP = 2;
%reduced_nGP = 1;
allNumel = [6; 12; 24; 48; 96];
timestepsize = 0.02;
totalTime = 1;

% data from Yu et al.
lengthBeam1 = sqrt(5);
lengthBeam2 = 1.5;
lengthBeam3 = sqrt(2);
E = 210000000000;
I1 = ((0.05)^4)/12;
I2 = I1;
I3 = I2;
A1 = 0.05^2;
A2 = A1;
A3 = A2;
rho1 = 2710 ;
rho2 = rho1 ;
rho3 = rho2 ;
nu = 0.3;
G = E / (2*(1+nu));

x = zeros(numel(allNumel),1);
y = x;

for index = 1:numel(allNumel)
    
    
    numberOfElements1 = allNumel(index)*10/6;
    numberOfElements2 = allNumel(index);
    numberOfElements3 = allNumel(index);
    timesteps = totalTime/timestepsize;
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = mfilename;
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = timesteps;
    setupObject.totalTime = totalTime;
    setupObject.newton.tolerance = 1e-8;
    setupObject.integrator = 'Midpoint';
    setupObject.plotObject.flag = true;
    setupObject.plotObject.steps = timesteps/5;
    setupObject.plotObject.keepFormerPlots = true;
    setupObject.plotObject.postPlotType = "stress";
    setupObject.plotObject.stress.component = 3; %select 3 for bending moment, 1 for normal force, and 2 for shear force
    setupObject.plotObject.view = [0,90];
    setupObject.plotObject.setXminXmax([-2.5;-2.5;-1],[4;2.5;1]);
    setupObject.plotObject.border = 0.1;
    setupObject.plotObject.plotInitialConfig = true;
    setupObject.plotObject.showGrid = true;
    setupObject.plotObject.colorscheme = 'viridis';
    
    dofObject = dofClass; % required object for dof and object handling
    dofObject.postDataObject.storeStateFlag = true;
    
    %% beam 1
    beamObject = beamClass(dofObject,2);
    beamObject.materialObject.name = "Hooke";
    beamObject.theory = beamType;
    beamObject.elementDisplacementType = 'mixedPH';
    %beamObject.elementDisplacementType = 'displacement';
    %beamObject.elementNameAdditionalSpecification = 'StanderStein';
    beamObject.materialObject.EI = E*I1;
    beamObject.materialObject.EA = E*A1;
    beamObject.materialObject.rhoI = rho1*I1;
    beamObject.materialObject.kGA = G*A1;
    beamObject.materialObject.rhoA = rho1*A1;
    beamObject.materialObject.rho = 1;
    
    % Spatial discretiaztion
    startpoint = [0;0];
    endpoint = [2;1];
    length = lengthBeam1;
    
    % Mesh
    [beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements1,order,startpoint,endpoint);
    beamObject.dimension = 1;
    beamObject.shapeFunctionObject.order = order;
    beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
    %beamObject.selectiveReducedShapeFunctionObject.order = order;
    %beamObject.selectiveReducedShapeFunctionObject.numberOfGausspoints = reduced_nGP;
    numberOfDispDOFFirstBeam = 3 * max(max(beamObject.meshObject.edof)); % displacements dofs
    numberOfDOFFirstBeam = numberOfDispDOFFirstBeam + numberOfElements1 * 3; % strain dofs
    %numberOfDOFFirstBeam = numberOfDispDOFFirstBeam;
    
    % dirichlet boundaries
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.nodeList = find(beamObject.meshObject.nodes(:, 1) == startpoint(1));
    dirichletObject1.nodalDof = [1, 2];
    dirichletObject1.masterObject = beamObject;
    
    %Bodyforce
    bodyForceObject = bodyForceClass(dofObject);
    bodyForceObject.typeOfLoad = 'deadLoad';
    bodyForceObject.masterObject = beamObject;
    bodyForceObject.loadFunction = rho1*A1*[0; -9.81];
    bodyForceObject.dimension = 2;
    bodyForceObject.timeFunction = str2func('@(t) 1');
    bodyForceObject.meshObject.edof = beamObject.meshObject.edof;
    bodyForceObject.shapeFunctionObject.order = 1;
    bodyForceObject.shapeFunctionObject.numberOfGausspoints = nGP;
    
    % solver
    beamObject.numericalTangentObject.computeNumericalTangent = true;
    beamObject.numericalTangentObject.showDifferences = false;
    
    %% beam 2
    beamObject2 = beamClass(dofObject,2);
    beamObject2.materialObject.name = "Hooke";
    beamObject2.theory = beamType;
    beamObject2.elementDisplacementType = 'mixedPH';
    %beamObject2.elementDisplacementType = 'displacement';
    %beamObject2.elementNameAdditionalSpecification = 'StanderStein';
    beamObject2.materialObject.EI = E*I2;
    beamObject2.materialObject.EA = E*A2;
    beamObject2.materialObject.rhoI = rho2*I2;
    beamObject2.materialObject.kGA = G*A2;
    beamObject2.materialObject.rhoA = rho2*A2;
    beamObject2.materialObject.rho = 1;
    
    % Spatial discretiaztion
    startpoint = [2;1];
    endpoint = [3.5;1];
    length = lengthBeam2;
    
    % Mesh
    [beamObject2.meshObject.nodes,beamObject2.meshObject.edof,edofNeumann2] = linearString(length,numberOfElements2,order,startpoint,endpoint);
    beamObject2.dimension = 1;
    beamObject2.shapeFunctionObject.order = order;
    beamObject2.shapeFunctionObject.numberOfGausspoints = nGP;
    %beamObject2.selectiveReducedShapeFunctionObject.order = order;
    %beamObject2.selectiveReducedShapeFunctionObject.numberOfGausspoints = reduced_nGP;
    
    numberOfDispDOFSecondBeam = 3 * max(max(beamObject2.meshObject.edof)); % displacements dofs
    numberOfDOFSecondBeam = numberOfDispDOFSecondBeam + numberOfElements2 * 3; % strain dofs
    %numberOfDOFSecondBeam = numberOfDispDOFSecondBeam;
    
    %Bodyforce
    bodyForceObject2 = bodyForceClass(dofObject);
    bodyForceObject2.typeOfLoad = 'deadLoad';
    bodyForceObject2.masterObject = beamObject2;
    bodyForceObject2.loadFunction = rho2*A2*[0; -9.81];
    bodyForceObject2.dimension = 2;
    bodyForceObject2.timeFunction = str2func('@(t) 1');
    bodyForceObject2.meshObject.edof = beamObject2.meshObject.edof;
    bodyForceObject2.shapeFunctionObject.order = 1;
    bodyForceObject2.shapeFunctionObject.numberOfGausspoints = nGP;
    
    % solver
    beamObject2.numericalTangentObject.computeNumericalTangent = true;
    beamObject2.numericalTangentObject.showDifferences = false;
    
    %% beam 3
    beamObject3 = beamClass(dofObject,2);
    beamObject3.materialObject.name = "Hooke";
    beamObject3.theory = beamType;
    beamObject3.elementDisplacementType = 'mixedPH';
    % beamObject3.elementDisplacementType = 'displacement';
    % beamObject3.elementNameAdditionalSpecification = 'StanderStein';
    beamObject3.materialObject.EI = E*I3;
    beamObject3.materialObject.EA = E*A3;
    beamObject3.materialObject.rhoI = rho3*I3;
    beamObject3.materialObject.kGA = G*A3;
    beamObject3.materialObject.rhoA = rho3*A3;
    beamObject3.materialObject.rho = 1;
    
    % Spatial discretiaztion
    startpoint = endpoint;
    endpoint = [2.5;0];
    length = lengthBeam3;
    
    % Mesh
    [beamObject3.meshObject.nodes,beamObject3.meshObject.edof,edofNeumann3] = linearString(length,numberOfElements3,order,startpoint,endpoint);
    beamObject3.dimension = 1;
    beamObject3.shapeFunctionObject.order = order;
    beamObject3.shapeFunctionObject.numberOfGausspoints = nGP;
    %beamObject3.selectiveReducedShapeFunctionObject.order = order;
    %beamObject3.selectiveReducedShapeFunctionObject.numberOfGausspoints = reduced_nGP;
    
    %Bodyforce
    bodyForceObject3 = bodyForceClass(dofObject);
    bodyForceObject3.typeOfLoad = 'deadLoad';
    bodyForceObject3.masterObject = beamObject3;
    bodyForceObject3.loadFunction = rho3*A3*[0; -9.81];
    bodyForceObject3.dimension = 2;
    bodyForceObject3.timeFunction = str2func('@(t) 1');
    bodyForceObject3.meshObject.edof = beamObject3.meshObject.edof;
    bodyForceObject3.shapeFunctionObject.order = 1;
    bodyForceObject3.shapeFunctionObject.numberOfGausspoints = nGP;
    
    % solver
    beamObject3.numericalTangentObject.computeNumericalTangent = true;
    beamObject3.numericalTangentObject.showDifferences = false;
    
    %dirichlet boundaries
    dirichletObject2 = dirichletClass(dofObject);
    dirichletObject2.nodeList = find(beamObject3.meshObject.nodes(:, 1) == endpoint(1));
    dirichletObject2.nodalDof = [1, 2];
    dirichletObject2.masterObject = beamObject3;
    
    %% Include constraints
    constraintObject1 = constraintClass(dofObject);
    constraintObject1.constrainedElement = numberOfElements1 + 1;
    constraintObject1.constrainedDOFinElement = [numberOfDOFFirstBeam + 1, numberOfDOFFirstBeam + 2];
    constraintObject1.constrainedDOFinOtherElement = [numberOfDispDOFFirstBeam - 2 , numberOfDispDOFFirstBeam - 1];
    
    constraintObject2 = constraintClass(dofObject);
    constraintObject2.constrainedElement = numberOfElements1 + numberOfElements2 + 1;
    constraintObject2.constrainedDOFinElement = [numberOfDOFFirstBeam + numberOfDOFSecondBeam + 1, numberOfDOFFirstBeam + numberOfDOFSecondBeam + 2];
    constraintObject2.constrainedDOFinOtherElement = [numberOfDOFFirstBeam + numberOfDispDOFSecondBeam - 2 , numberOfDOFFirstBeam + numberOfDispDOFSecondBeam - 1];
    
    %% Solve!
    dofObject = runNewton(setupObject, dofObject);
    
    x(index) = dofObject.postDataObject.stateJournal(end).position(numberOfDOFFirstBeam+3*(numberOfElements2/2)+1);
    y(index) = dofObject.postDataObject.stateJournal(end).position(numberOfDOFFirstBeam+3*(numberOfElements2/2)+2);
    
end

x_ref = 1.399354522584749;
y_ref = -0.757502351760989;

error = zeros(numel(allNumel)-1,1);
for index2 = 1:numel(allNumel)-1
    error(index2) = norm([x(index2);y(index2)]-[x_ref;y_ref])/norm([x_ref;y_ref]);
end
figure();
loglog(lengthBeam2./allNumel(1:end-1),error);
%matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'space_convergence.tikz', 'showInfo', false, 'floatformat', '%.4g')
