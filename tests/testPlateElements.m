% unit test for plate elements

disp('=======================================');
disp('Now testing: plate elements');
disp('=======================================');

elementsToTest = {'displacementHookeEndpoint', ...
    'displacementHookeEndpointQ8', ...
    'displacementHookeEndpointQ9', ...
    'displacementHughesTezduyarHookeEndpoint'};

for ii = 1:size(elementsToTest, 2)
    displacementType = 'displacement';
    materialName = 'Hooke';
    order = 1;
    serendipity = false;
    elementNameAdditionalSpecification = '';
    switch elementsToTest{ii}
        case 'displacementHookeEndpoint'
            correctSolution = 0.37301;
        case 'displacementHookeEndpointQ8'
            order = 2;
            serendipity = true;
            correctSolution = 26.9312;
        case 'displacementHookeEndpointQ9'
            order = 2;
            correctSolution = 57.2402;
        case 'displacementHughesTezduyarHookeEndpoint'
            elementNameAdditionalSpecification = 'HughesTezduyar';
            correctSolution = 57.0195;
        case 'displacementBatheDvorkinHookeEndpoint'
            elementNameAdditionalSpecification = 'BatheDvorkin';
            correctSolution = 56.9754;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    disp('---------------------------------------');
    disp(['- Current Element: ', elementsToTest{ii}]);
    disp('---------------------------------------');

    % setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'plateTest';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';

    dofObject = dofClass; % required object for dof and object handling

    % continuum Objects
    plateLengthX = 50;
    plateLengthY = 50;
    plateObject = plateClass(dofObject);
    plateObject.elementDisplacementType = displacementType;
    plateObject.elementNameAdditionalSpecification = elementNameAdditionalSpecification;
    [plateObject.meshObject.nodes, plateObject.meshObject.edof] = meshPatchTestDistorted2D(plateLengthX, plateLengthY, order, serendipity);
    plateObject.materialObject.name = materialName;
    plateObject.materialObject.rho = 0;
    plateObject.materialObject.E = 1e4;
    plateObject.materialObject.nu = 0.3;
    plateObject.h = 1;
    plateObject.dimension = 2;
    plateObject.shapeFunctionObject.order = order;
    plateObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

    % Dirichlet boundary
    boundary1 = dirichletClass(dofObject);
    boundary1.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
    boundary1.nodalDof = 1:3;
    boundary1.masterObject = plateObject;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.masterObject = plateObject;
    neumannObject.loadGeometry = 'area';
    neumannObject.loadVector = -[1e-2; 0; 0];
    neumannObject.meshObject.edof = plateObject.meshObject.edof;

    % Nodal Forces
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = plateObject;
    nodalLoadObject.loadVector = -[50; 0; 0];
    nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:, 1) == plateLengthX & plateObject.meshObject.nodes(:, 2) == plateLengthY);

    % runNewton
    dofObject = runNewton(setupObject, dofObject);
    maxDisplacement = max(abs(plateObject.qN1(:, 1)));
    disp(['Maximum displacement: ', num2str(maxDisplacement)]);
    assert(abs(maxDisplacement-correctSolution) < 10^-4, ['Wrong results for ', elementsToTest{ii}]);
end