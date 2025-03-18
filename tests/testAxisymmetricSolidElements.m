% unit test for axisymmetric solid elements

disp('=======================================');
disp('Now testing: axisymmetric solid elements');
disp('=======================================');

elementsToTest = {'displacementHookeEndpoint', ...
    'displacementHookeEndpointQ9', ...
    'displacementPetrovGalerkinHookeEndpoint'};

for ii = 1:size(elementsToTest, 2)
    displacementType = 'displacement';
    materialName = 'Hooke';
    order = 1;
    serendipity = false;
    elementNameAdditionalSpecification = '';
    mixedTypeShapeFunctionData = 0;
    condensation = false;
    switch elementsToTest{ii}
        case 'displacementHookeEndpoint'
            correctSolution = 0.0068352;
        case 'displacementHookeEndpointQ9'
            order = 2;
            correctSolution = 0.0070348;
        case 'displacementPetrovGalerkinHookeEndpoint'
            elementNameAdditionalSpecification = 'PetrovGalerkin';
            correctSolution = 0.0068335;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    disp('---------------------------------------');
    disp(['- Current Element: ', elementsToTest{ii}]);
    disp('---------------------------------------');

    % setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'axisymmetricSolidTest';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';

    dofObject = dofClass; % required object for dof and object handling

    % continuum objects
    axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
    axisymmetricSolidObject.dimension = 2;
    
    lengthX = 1;
    lengthY = 1;
    numberOfElementsX = 2;
    numberOfElementsY = 2;
    innerRadius = 1;
    axisymmetricSolidObject.shapeFunctionObject.order = order;
    axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
    
    [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
    nodeToDistort = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
    axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;
    axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) + 0.12;
    axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) + 0.11;
    axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);

    axisymmetricSolidObject.elementDisplacementType = displacementType;
    axisymmetricSolidObject.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
    
    axisymmetricSolidObject.materialObject.name = materialName;
    axisymmetricSolidObject.materialObject.rho = 0;
    axisymmetricSolidObject.materialObject.E = 250;
    axisymmetricSolidObject.materialObject.nu = 0.1;

    axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = mixedTypeShapeFunctionData;
    axisymmetricSolidObject.mixedFEObject.condensation = condensation;
    
    % Dirichlet boundary
    dirichletBoundary = dirichletClass(dofObject);
    dirichletBoundary.nodeList = find(abs(axisymmetricSolidObject.meshObject.nodes(:, 2) + lengthY/2) < 1e-10);
    dirichletBoundary.nodalDof = 2;
    dirichletBoundary.timeFunction = @(t, XYZ) -lengthY/2;
    dirichletBoundary.masterObject = axisymmetricSolidObject;
    
    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(abs(axisymmetricSolidObject.meshObject.nodes(:, 2) - lengthY/2) < 1e-10);
    dirichletBoundary2.nodalDof = 2;
    dirichletBoundary2.timeFunction = @(t, XYZ) lengthY/2;
    dirichletBoundary2.masterObject = axisymmetricSolidObject;
    
    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = axisymmetricSolidObject;
    neumannObject.loadVector = [1; 0];
    neumannObject.meshObject.edof = bounEdof.SX1;

    % runNewton
    dofObject = runNewton(setupObject, dofObject);
    nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
    displacement = axisymmetricSolidObject.qN1(nodeToMeasure, 1) - axisymmetricSolidObject.qN(nodeToMeasure, 1);
    disp(['Computed displacement: ', num2str(displacement)]);
    assert(abs(displacement-correctSolution) < 10^-7, ['Wrong results for ', elementsToTest{ii}]);
end