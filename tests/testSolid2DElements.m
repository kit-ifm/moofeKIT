% unit test for 2D solid elements (static test)

disp('=======================================');
disp('Now testing: 2D solid elements');
disp('=======================================');

elementsToTest = {'displacementHooke2DEndpoint', ...
    'displacementHooke2DEndpointQ8', ...
    'displacementHooke2DEndpointQ9', ...
    'displacementNeoHooke2DEndpoint', ...
    'displacementSCSaintVenant2DEndpoint', ...
    'displacementPetrovGalerkinRajendranLiew2003Hooke2DEndpoint', ...
    'displacementPetrovGalerkinXieEtAl2016UQ8Hooke2DEndpoint', ...
    'displacementPetrovGalerkinCenEtAl2015Hooke2DEndpoint', ...
    'displacementPetrovGalerkinXieEtAl2016TQ4Hooke2DEndpoint', ...
    'displacementPetrovGalerkinSaintVenant2DEndpoint', ...
    'selectiveReducedIntegrationHooke2DEndpoint', ...
    'incompatibleModesPetrovGalerkinHuangEtAl2020Hooke2DEndpoint', ...
    'pianSumiharaHooke2DEndpoint', ...
    'pianSumiharaSaintVenant2DEndpoint', ...
    'easHooke2DEndpoint', ...
    'easPetrovGalerkinPfefferkornBetsch2021Hooke2DEndpoint'};

for ii = 1:size(elementsToTest, 2)
    elementProperties = getElementProperties(elementsToTest{ii});

    disp('---------------------------------------');
    disp(['- Current Element: ', elementsToTest{ii}]);
    disp('---------------------------------------');

    % setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = '2DSolidTest';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = elementProperties.numberOfLoadSteps;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.newton.tolerance = 1e-4;
    setupObject.integrator = elementProperties.integrator;

    dofObject = dofClass; % required object for dof and object handling

    % continuum Objects
    solidObject = solidClass(dofObject);
    solidObject.elementDisplacementType = elementProperties.displacementType;
    solidObject.elementNameAdditionalSpecification = elementProperties.elementNameAdditionalSpecification;
    [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshCooksMembrane(4, 4, elementProperties.order, elementProperties.serendipity);
    solidObject.materialObject.name = elementProperties.materialName;
    solidObject.elementNameAdditionalSpecification = elementProperties.elementNameAdditionalSpecification;
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.E = 1;
    solidObject.materialObject.nu = 1/3;
    E = solidObject.materialObject.E;
    nu = solidObject.materialObject.nu;
    solidObject.materialObject.lambda = nu / (1 - 2 * nu) * (1 / (1 + nu)) * E;
    solidObject.materialObject.mu = 1 / 2 * (1 / (1 + nu)) * E;
    solidObject.dimension = 2;
    solidObject.shapeFunctionObject.order = elementProperties.order;
    solidObject.shapeFunctionObject.numberOfGausspoints = elementProperties.numberOfGausspoints(elementProperties.order);
    solidObject.mixedFEObject.typeShapeFunctionData = elementProperties.mixedFEObjectTypeShapeFunctionData;
    solidObject.mixedFEObject.condensation = true;

    % Dirichlet boundary
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.masterObject = solidObject;
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
    dirichletObject.nodalDof = 1:2;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = solidObject;
    neumannObject.loadVector = [0; 1/16];
    neumannObject.timeFunction = @(t) t;
    neumannObject.meshObject.edof = edofNeumann;

    % runNewton
    dofObject = runNewton(setupObject, dofObject);
    nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-48) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-52) < 1e-8); % middle right node
    yDisplacement = solidObject.qN1(nodeToMeasure,2)-solidObject.qR(nodeToMeasure,2);
    disp(['Displacement: ', num2str(yDisplacement)]);
    assert(abs(yDisplacement-elementProperties.correctSolution) < 10^-4, ['Wrong results for ', elementsToTest{ii}]);
end

function elementProperties = getElementProperties(elementName)
displacementType = 'displacement';
materialName = 'HookeESZ';
order = 1;
serendipity = false;
numberOfGausspoints = @(order) (order + 1)^2;
elementNameAdditionalSpecification = '';
mixedFEObjectTypeShapeFunctionData = 0;
numberOfLoadSteps = 1;
integrator = 'Endpoint';
switch elementName
    case 'displacementHooke2DEndpoint'
        correctSolution = 18.2992;
    case 'displacementHooke2DEndpointQ8'
        order = 2;
        serendipity = true;
        correctSolution = 23.7083;
    case 'displacementHooke2DEndpointQ9'
        order = 2;
        correctSolution = 23.8397;
    case 'displacementNeoHooke2DEndpoint'
        materialName = 'NeoHookeEVZ'; % TODO: change to ESZ
        numberOfLoadSteps = 5;
        correctSolution = 12.6413;
    case 'displacementSCSaintVenant2DEndpoint'
        elementNameAdditionalSpecification = 'SC';
        materialName = 'SaintVenantESZ';
        numberOfLoadSteps = 5;
        correctSolution = 13.5356;
    case 'displacementPetrovGalerkinRajendranLiew2003Hooke2DEndpoint'
        order = 2;
        serendipity = true;
        elementNameAdditionalSpecification = 'PetrovGalerkinRajendranLiew2003';
        correctSolution = 23.7535;
    case 'displacementPetrovGalerkinCenEtAl2015Hooke2DEndpoint'
        elementNameAdditionalSpecification = 'PetrovGalerkinCenEtAl2015';
        correctSolution = 23.4337;
    case 'displacementPetrovGalerkinXieEtAl2016UQ8Hooke2DEndpoint'
        order = 2;
        serendipity = true;
        elementNameAdditionalSpecification = 'PetrovGalerkinXieEtAl2016UQ8';
        correctSolution = 23.7257;
    case 'displacementPetrovGalerkinXieEtAl2016TQ4Hooke2DEndpoint'
        elementNameAdditionalSpecification = 'PetrovGalerkinXieEtAl2016TQ4';
        correctSolution = 23.4337;
    case 'displacementPetrovGalerkinSaintVenant2DEndpoint'
        elementNameAdditionalSpecification = 'PetrovGalerkin';
        materialName = 'SaintVenantESZ';
        numberOfLoadSteps = 5;
        correctSolution = 13.6742;
    case 'selectiveReducedIntegrationHooke2DEndpoint'
        displacementType = 'selectiveReducedIntegration';
        materialName = 'HookeEVZ'; % TODO: change to ESZ
        correctSolution = 17.7377;
    case 'incompatibleModesPetrovGalerkinHuangEtAl2020Hooke2DEndpoint'
        displacementType = 'incompatibleModes';
        elementNameAdditionalSpecification = 'PetrovGalerkinHuangEtAl2020';
        mixedFEObjectTypeShapeFunctionData = 10;
        numberOfGausspoints = @(order) (order + 2)^2;
        correctSolution = 23.2899;
    case 'pianSumiharaHooke2DEndpoint'
        displacementType = 'pianSumihara';
        mixedFEObjectTypeShapeFunctionData = 5;
        correctSolution = 23.0219;
    case 'pianSumiharaSaintVenant2DEndpoint'
        displacementType = 'pianSumihara';
        materialName = 'SaintVenantESZ';
        mixedFEObjectTypeShapeFunctionData = 5;
        numberOfLoadSteps = 5;
        correctSolution = 14.0669;
    case 'easHooke2DEndpoint'
        displacementType = 'eas';
        mixedFEObjectTypeShapeFunctionData = 4;
        correctSolution = 23.0164;
    case 'easPetrovGalerkinPfefferkornBetsch2021Hooke2DEndpoint'
        displacementType = 'eas';
        elementNameAdditionalSpecification = 'PetrovGalerkinPfefferkornBetsch2021';
        mixedFEObjectTypeShapeFunctionData = 4;
        correctSolution = 23.4337;
    otherwise
        warning(['Element ', elementName, ' not implemented!'])
end

elementProperties = struct();
elementProperties.displacementType = displacementType;
elementProperties.materialName = materialName;
elementProperties.order = order;
elementProperties.serendipity = serendipity;
elementProperties.numberOfGausspoints = numberOfGausspoints;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
elementProperties.mixedFEObjectTypeShapeFunctionData = mixedFEObjectTypeShapeFunctionData;
elementProperties.numberOfLoadSteps = numberOfLoadSteps;
elementProperties.integrator = integrator;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;

elementProperties.correctSolution = correctSolution;
end