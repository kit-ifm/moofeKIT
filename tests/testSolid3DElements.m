% unit test for 3D solid elements (static test)

disp('=======================================');
disp('Now testing: 3D solid elements');
disp('=======================================');

elementsToTest = {'displacementHookeEndpoint', ...
    'displacementHookeEndpointH20', ...
    'displacementHookeEndpointH27', ...
    'displacementHookeMidpoint', ...
    'displacementHookeSplitMidpoint', ...
    'displacementSCSaintVenantEndpoint', ...
    'displacementSCSaintVenantMidpoint', ...
    'easHookeEndpoint', ...
    'easPetrovGalerkinPfefferkornBetsch2021HookeEndpoint', ...
    'displacementPetrovGalerkinZhouEtAl2017HookeEndpoint', ...
    'displacementPetrovGalerkinXieEtAl2016UH20HookeEndpoint', ...
    'displacementSCMooneyRivlinEndpoint', ...
    'displacementSCMooneyRivlinMidpoint', ...
    'displacementSCMooneyRivlinDiscreteGradient', ...
    'mixedSCMooneyRivlinEndpointCascadeCGc', ...
    'mixedSCMooneyRivlinMidpoint', ...
    'mixedSCMooneyRivlinDiscreteGradient'};
%     'mixedSCMooneyRivlinEndpointNonCascadeCGc', ... % FIXME: analytical tangent not implemented
%     'mixedSCMooneyRivlinEndpointInvariantCGc', ...  % FIXME: analytical tangent not implemented


for ii = 1:size(elementsToTest, 2)
    elementProperties = getElementProperties(elementsToTest{ii});

    disp('---------------------------------------');
    disp(['- Current Element: ', elementsToTest{ii}]);
    disp('---------------------------------------');

    % setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = '3DSolidTest';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = elementProperties.numberOfLoadSteps;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.newton.tolerance = 1e-4;
    setupObject.integrator = elementProperties.integrator;

    dofObject = dofClass; % required object for dof and object handling

    % continuum Objects
    nel = 4;
    solidObject = solidClass(dofObject);
    solidObject.elementDisplacementType = elementProperties.displacementType;
    solidObject.elementNameAdditionalSpecification = elementProperties.elementNameAdditionalSpecification;
    [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshGeneratorCooksMembrane(nel, 1, nel, elementProperties.order, elementProperties.serendipity);
    solidObject.materialObject.name = elementProperties.materialName;
    solidObject.elementNameAdditionalSpecification = elementProperties.elementNameAdditionalSpecification;
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.E = 2261;
    solidObject.materialObject.nu = 0.4955;
    E = solidObject.materialObject.E;
    nu = solidObject.materialObject.nu;
    solidObject.materialObject.lambda = nu / (1 - 2 * nu) * (1 / (1 + nu)) * E;
    solidObject.materialObject.mu = 1 / 2 * (1 / (1 + nu)) * E;
    solidObject.materialObject.a = solidObject.materialObject.mu / 6;
    solidObject.materialObject.b = solidObject.materialObject.mu / 3;
    solidObject.materialObject.c = solidObject.materialObject.lambda;
    solidObject.materialObject.d = 2 * (solidObject.materialObject.a + 2 * solidObject.materialObject.b);
    solidObject.dimension = 3;
    solidObject.shapeFunctionObject.order = elementProperties.order;
    solidObject.shapeFunctionObject.numberOfGausspoints = (elementProperties.order + 1)^3;
    solidObject.mixedFEObject.typeShapeFunctionData = elementProperties.mixedFEObjectTypeShapeFunctionData;
    solidObject.mixedFEObject.condensation = true;

    % Dirichlet boundary
    boundary1 = dirichletClass(dofObject);
    boundary1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
    boundary1.nodalDof = 1:3;
    boundary1.masterObject = solidObject;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.masterObject = solidObject;
    neumannObject.loadVector = [0; 0; 100];
    neumannObject.timeFunction = @(t) t;
    neumannObject.loadGeometry = 'area';
    neumannObject.meshObject.edof = edofBoundary.SX2;

    % runNewton
    dofObject = runNewton(setupObject, dofObject);
    maxDisplacement = max(abs(solidObject.qN1(end, 3)-solidObject.qR(end, 3)));
    disp(['Maximum displacement: ', num2str(maxDisplacement)]);
    assert(abs(maxDisplacement-elementProperties.correctSolution) < 10^-4, ['Wrong results for ', elementsToTest{ii}]);
end

function elementProperties = getElementProperties(elementName)
displacementType = 'displacement';
materialName = 'Hooke';
order = 1;
serendipity = false;
elementNameAdditionalSpecification = '';
mixedFEObjectTypeShapeFunctionData = 0;
numberOfLoadSteps = 1;
integrator = 'Endpoint';
switch elementName
    case 'displacementHookeEndpoint'
        correctSolution = 5.8617;
    case 'displacementHookeEndpointH20'
        order = 2;
        serendipity = true;
        correctSolution = 16.198;
    case 'displacementHookeEndpointH27'
        order = 2;
        correctSolution = 16.3295;
    case 'displacementHookeMidpoint'
        integrator = 'Midpoint';
        correctSolution = 5.8617;
    case 'displacementHookeSplitMidpoint'
        integrator = 'Midpoint';
        correctSolution = 5.8617;
    case 'displacementSCSaintVenantEndpoint'
        correctSolution = 5.8617;
    case 'displacementSCSaintVenantMidpoint'
        integrator = 'Midpoint';
        correctSolution = 5.8617;
    case 'easHookeEndpoint'
        displacementType = 'eas';
        mixedFEObjectTypeShapeFunctionData = 9;
        correctSolution = 15.4177;
    case 'easPetrovGalerkinPfefferkornBetsch2021HookeEndpoint'
        displacementType = 'eas';
        elementNameAdditionalSpecification = 'PetrovGalerkinPfefferkornBetsch2021';
        mixedFEObjectTypeShapeFunctionData = 12;
        correctSolution = 16.1567;
    case 'displacementPetrovGalerkinZhouEtAl2017HookeEndpoint'
        elementNameAdditionalSpecification = 'PetrovGalerkinZhouEtAl2017';
        correctSolution = 15.6864;
    case 'displacementPetrovGalerkinXieEtAl2016UH20HookeEndpoint'
        order = 2;
        serendipity = true;
        elementNameAdditionalSpecification = 'PetrovGalerkinXieEtAl2016UH20';
        correctSolution = 16.2306;
    case 'displacementSCMooneyRivlinEndpoint'
        displacementType = 'displacementSC';
        materialName = 'MooneyRivlin';
        numberOfLoadSteps = 8;
        order = 2;
        serendipity = true;
        correctSolution = 11.8012;
    case 'displacementSCMooneyRivlinMidpoint'
        displacementType = 'displacementSC';
        materialName = 'MooneyRivlin';
        order = 2;
        numberOfLoadSteps = 8;
        integrator = 'Midpoint';
        serendipity = true;
        correctSolution = 11.7916;
    case 'displacementSCMooneyRivlinDiscreteGradient'
        displacementType = 'displacementSC';
        materialName = 'MooneyRivlin';
        order = 2;
        numberOfLoadSteps = 8;
        integrator = 'DiscreteGradient';
        serendipity = true;
        correctSolution = 11.8035;
    case 'mixedSCMooneyRivlinEndpointCascadeCGc'
        displacementType = 'mixedSC';
        materialName = 'MooneyRivlin';
        order = 2;
        serendipity = true;
        mixedFEObjectTypeShapeFunctionData = 1;
        correctSolution = 12.3679;
    case 'mixedSCMooneyRivlinEndpointNonCascadeCGc'
        displacementType = 'mixedSC';
        materialName = 'MooneyRivlin';
        order = 2;
        serendipity = true;
        mixedFEObjectTypeShapeFunctionData = 1;
        elementNameAdditionalSpecification = 'NonCascadeCGc';

        correctSolution = 12.3679;
    case 'mixedSCMooneyRivlinEndpointInvariantCGc'
        displacementType = 'mixedSC';
        materialName = 'MooneyRivlin';
        order = 2;
        serendipity = true;
        mixedFEObjectTypeShapeFunctionData = 1;
        elementNameAdditionalSpecification = 'InvariantCGc';

        correctSolution = 12.3679;
    case 'mixedSCMooneyRivlinMidpoint'
        displacementType = 'mixedSC';
        materialName = 'MooneyRivlin';
        order = 2;
        serendipity = true;
        mixedFEObjectTypeShapeFunctionData = 1;
        integrator = 'Midpoint';

        correctSolution = 14.5902;
    case 'mixedSCMooneyRivlinDiscreteGradient'
        displacementType = 'mixedSC';
        materialName = 'MooneyRivlin';
        order = 2;
        serendipity = true;
        mixedFEObjectTypeShapeFunctionData = 1;
        integrator = 'DiscreteGradient';

        correctSolution = 12.9262;
    otherwise
        warning(['Element ', elementName, ' not implemented!'])
end

elementProperties = struct();
elementProperties.displacementType = displacementType;
elementProperties.materialName = materialName;
elementProperties.order = order;
elementProperties.serendipity = serendipity;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
elementProperties.mixedFEObjectTypeShapeFunctionData = mixedFEObjectTypeShapeFunctionData;
elementProperties.numberOfLoadSteps = numberOfLoadSteps;
elementProperties.integrator = integrator;
elementProperties.elementNameAdditionalSpecification = elementNameAdditionalSpecification;

elementProperties.correctSolution = correctSolution;
end