% unit test for geometrically exact beam elements 

disp('=======================================');
disp('Now testing: Geometrically exact beam elements');
disp('=======================================');

elementsToTest = {'displacementStanderSteinGeometricallyExactHookeMidpoint', ...
    'displacementPhIrreducibleGeometricallyExactHookeMidpoint', ...
    'mixedPHGeometricallyExactHookeEndpoint', ...
    'mixedPHGeometricallyExactHookeMidpoint'};

for ii = 1:size(elementsToTest, 2)
    elementProperties = getElementProperties(elementsToTest{ii});

    disp('---------------------------------------');
    disp(['- Current Element: ', elementsToTest{ii}]);
    disp('---------------------------------------');
    
    rotationAngle = 0;
    
    %% Parameters
    lengthBeam = 1;
    beamType = 'GeometricallyExact';
    numberOfElements = 10;
    E = 184000;
    R = 0.186;
    selRednGP = 1;
    rho = 920;
    A = pi*R^2;
    I = (pi*R^4)/4;
    rhoA = rho*A;
    shearCorFact = 6/7; %shear correction factor (rectangle 5/6, circle 6/7)
    
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = mfilename;
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 2;
    setupObject.totalTime = 0.2;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.newton.tolerance = 1e-10;
    setupObject.integrator = elementProperties.integrator;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    
    dofObject = dofClass; % required object for dof and object handling
    dofObject.postDataObject.storeStateFlag = true;
    
    %% continuum Objects
    beamObject = beamClass(dofObject,2);
    beamObject.materialObject.name = "Hooke";
    beamObject.theory = beamType;
    beamObject.elementDisplacementType = elementProperties.displacementType;
    beamObject.elementNameAdditionalSpecification = elementProperties.elementNameAdditionalSpecification;
    beamObject.materialObject.E = E;
    nu = 0.3;
    beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
    beamObject.materialObject.I = I;
    beamObject.materialObject.A = A;
    beamObject.materialObject.rho = rho;
    beamObject.materialObject.shearCorrectionCoefficient = shearCorFact;
    
    startpoint = [0;0];
    endpoint = [cos(rotationAngle)*lengthBeam;sin(rotationAngle)*lengthBeam];
    length = lengthBeam;
    
    % Mesh
    [beamObject.meshObject.nodes,beamObject.meshObject.edof,edofNeumann] = linearString(length,numberOfElements,1,startpoint,endpoint);
    beamObject.dimension = 1;
    beamObject.shapeFunctionObject.order = 1;
    beamObject.shapeFunctionObject.numberOfGausspoints = 2;
    beamObject.selectiveReducedShapeFunctionObject.order = 1;
    beamObject.selectiveReducedShapeFunctionObject.numberOfGausspoints = selRednGP;
    
    % dirichlet boundaries
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = [1, 2, 3];
    dirichletObject.masterObject = beamObject;
    
    % neumann boundaries
    Q = 1000*[-sin(rotationAngle);cos(rotationAngle);0];
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = beamObject;
    nodalLoadObject.loadVector = Q;
    nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == endpoint(1));
    nodalLoadObject.timeFunction = str2func('@(t) sin(pi*t/0.2).*(t<=0.2+1e-12)');

    %Bodyforce
    bodyForceObject = bodyForceClass(dofObject);
    bodyForceObject.typeOfLoad = 'deadLoad';
    bodyForceObject.masterObject = beamObject;
    bodyForceObject.loadFunction = rhoA*[0; -9.81];
    bodyForceObject.dimension = 2;
    bodyForceObject.timeFunction = str2func('@(t) 1');
    bodyForceObject.meshObject.edof = beamObject.meshObject.edof;
    bodyForceObject.shapeFunctionObject.order = 1;
    bodyForceObject.shapeFunctionObject.numberOfGausspoints = 2;
    
    %% solver
    beamObject.numericalTangentObject.computeNumericalTangent = true;
    beamObject.numericalTangentObject.showDifferences = false;
    
    % runNewton
    dofObject = runNewton(setupObject, dofObject);
    nodeToMeasure = numberOfElements+1;
    DisplacementAndRotation = beamObject.qN1(nodeToMeasure,:)-beamObject.qR(nodeToMeasure,:);
    disp(['Displacement and Rotation: ', num2str(DisplacementAndRotation)]);
    assert(~any((DisplacementAndRotation-elementProperties.correctSolution) > 10^-5*[1,1,1]) , ['Wrong results for ', elementsToTest{ii}]);
end

function elementProperties = getElementProperties(elementName)
displacementType = 'displacement';
elementNameAdditionalSpecification = '';
integrator = 'Midpoint';
switch elementName
    case 'displacementStanderSteinGeometricallyExactHookeMidpoint'
        elementNameAdditionalSpecification = 'StanderStein';
        correctSolution = [-0.093072 0.29187 0.74402];
    case 'displacementPhIrreducibleGeometricallyExactHookeMidpoint'
        elementNameAdditionalSpecification = 'PhIrreducible';
        correctSolution = [-0.08888 0.28768 0.73163];
    case 'mixedPHGeometricallyExactHookeEndpoint'
        displacementType = 'mixedPH';
        integrator = 'Endpoint';
        correctSolution = [-0.07859 0.16118 0.33823];
    case 'mixedPHGeometricallyExactHookeMidpoint'
        displacementType = 'mixedPH';
        correctSolution = [-0.091077 0.28976 0.74456];
    otherwise
        warning(['Element ', elementName, ' not implemented!'])
end

elementProperties = struct();
elementProperties.displacementType = displacementType;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
elementProperties.integrator = integrator;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
elementProperties.correctSolution = correctSolution;

end