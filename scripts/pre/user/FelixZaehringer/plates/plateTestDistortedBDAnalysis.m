%% Start parameters
elementsToTest = {'displacementPetrovGalerkinBatheDvorkin', 'displacementBatheDvorkin'}; %, 'dispQ8', 'dispPetrovGalerkin', 'BD', 'easBD', 'eas'}; % {'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD'};

%% Script
resultArray = zeros(length(elementsToTest));
s = 1;

for ii = 1:length(elementsToTest)

    disp('=======================================');
    disp(['Current element: ', elementsToTest{ii}]);
    disp(['Current s: ', num2str(s)]);
    disp('=======================================');

    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'plateSimple';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.usePreconditioning = false;

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    plateObject = plateClass(dofObject);
    switch elementsToTest{ii}
        case 'displacementQ8'
            plateObject.materialObject.name = 'Hooke';
            order = 2;
            serendipity = true;
        case 'displacementQ4'
            plateObject.materialObject.name = 'Hooke';
            order = 1;
            serendipity = false;
        case 'selectiveReducedIntegration'
            plateObject.materialObject.name = 'Hooke';
            order = 1;
            plateObject.elementDisplacementType = 'selectiveReducedIntegration';
            serendipity = false;
            plateObject.mixedFEObject.condensation = true;
        case 'selectiveReducedIntegrationPetrovGalerkin'
            plateObject.materialObject.name = 'Hooke';
            order = 1;
            plateObject.elementDisplacementType = 'selectiveReducedIntegration';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
            serendipity = false;
            plateObject.mixedFEObject.condensation = true;
        case 'displacementPetrovGalerkinQ4'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
            order = 1;
            serendipity = false;
        case 'displacementQ9'
            plateObject.materialObject.name = 'Hooke';
            order = 2;
            serendipity = false;
        case 'displacementQ16'
            plateObject.materialObject.name = 'Hooke';
            order = 3;
            serendipity = false;
        case 'displacementPetrovGalerkinQ8'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementDisplacementType = 'displacement';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
            order = 2;
            serendipity = true;
        case 'displacementPetrovGalerkinBatheDvorkin'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementDisplacementType = 'displacement';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin';
            order = 1;
            serendipity = false;
        case 'displacementBatheDvorkin'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
            order = 1;
            serendipity = false;
        case 'displacementHughesTezduyar'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementNameAdditionalSpecification = 'HughesTezduyar';
            order = 1;
            serendipity = false;
        case 'displacementPetrovGalerkinHughesTezduyar'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinHughesTezduyar';
            order = 1;
            serendipity = false;
        case 'easSimoRifai'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementDisplacementType = 'eas';
            plateObject.elementNameAdditionalSpecification = 'SimoRifai';
            plateObject.mixedFEObject.condensation = true;
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            order = 1;
            serendipity = false;
        case 'easPetrovGalerkinSimoRifai'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementDisplacementType = 'easPetrovGalerkin';
            plateObject.elementNameAdditionalSpecification = 'SimoRifai';
            plateObject.mixedFEObject.condensation = true;
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            order = 1;
            serendipity = false;
        case 'easAndelfingerRamm'
            plateObject.materialObject.name = 'HookeBD';
            plateObject.elementDisplacementType = 'eas';
            plateObject.mixedFEObject.condensation = true;
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            order = 1;
            serendipity = false;
        case 'easPetrovGalerkin'
            plateObject.materialObject.name = 'Hooke';
            plateObject.elementDisplacementType = 'easPetrovGalerkin';
            plateObject.mixedFEObject.condensation = true;
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            order = 1;
            serendipity = false;
        case 'easPetrovGalerkinAndelfingerRamm'
            plateObject.materialObject.name = 'HookeBD';
            plateObject.elementDisplacementType = 'easPetrovGalerkin';
            plateObject.mixedFEObject.condensation = true;
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            order = 1;
            serendipity = false;
        case 'easPetrovGalerkinSelectiveReducedIntegration'
            plateObject.materialObject.name = 'Hooke';
            order = 1;
            plateObject.elementDisplacementType = 'eas';
            plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinSelectiveReducedIntegration';
            plateObject.mixedFEObject.typeShapeFunctionData = 4;
            serendipity = false;
            plateObject.mixedFEObject.condensation = true;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    plateLengthX = 2;
    plateLengthY = 4;

    [plateObject.meshObject.nodes, plateObject.meshObject.edof, bounEdof] = meshRectangle(plateLengthX, plateLengthY, 1, 2, order, serendipity);
    nodeToDistort1 = find(plateObject.meshObject.nodes(:, 1) == -plateLengthX/2 & plateObject.meshObject.nodes(:, 2) == 0);
    nodeToDistort2 = find(plateObject.meshObject.nodes(:, 1) == plateLengthX/2 & plateObject.meshObject.nodes(:, 2) == 0);
    plateObject.meshObject.nodes(nodeToDistort1, 2) = plateObject.meshObject.nodes(nodeToDistort1, 2) + s;
    plateObject.meshObject.nodes(nodeToDistort2, 2) = plateObject.meshObject.nodes(nodeToDistort2, 2) - s;
    if order == 2 && serendipity
        plateObject.meshObject.nodes(nodeToDistort1-2, 2) = plateObject.meshObject.nodes(nodeToDistort1-2, 2) + s / 2;
        plateObject.meshObject.nodes(nodeToDistort1+3, 2) = plateObject.meshObject.nodes(nodeToDistort1+3, 2) + s / 2;
        plateObject.meshObject.nodes(nodeToDistort2-3, 2) = plateObject.meshObject.nodes(nodeToDistort2-3, 2) - s / 2;
        plateObject.meshObject.nodes(nodeToDistort2+2, 2) = plateObject.meshObject.nodes(nodeToDistort2+2, 2) - s / 2;
    end
    plateObject.meshObject.nodes(:, 1) = plateObject.meshObject.nodes(:, 1) + plateLengthX / 2;
    plateObject.meshObject.nodes(:, 2) = plateObject.meshObject.nodes(:, 2) + plateLengthY / 2;

    %         [plateObject.meshObject.nodes,plateObject.meshObject.edof, bounEdof] = meshRectangle(plateLengthX, plateLengthY, 80, 80, order, serendipity);
    %         plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
    %         plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;

    %plateObject.materialObject.name = 'Hooke';
    plateObject.materialObject.rho = 0;
    plateObject.materialObject.E = 6;
    plateObject.materialObject.nu = 0;
    plateObject.h = 1;
    plateObject.dimension = 2;
    plateObject.shapeFunctionObject.order = order;
    plateObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

    % Dirichlet boundary
    boundary1 = dirichletClass(dofObject);
    boundary1.nodeList = find(plateObject.meshObject.nodes(:, 2) == plateLengthY);
    boundary1.nodalDof = 1:3;
    boundary1.masterObject = plateObject;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = plateObject;
    neumannObject.loadVector = 1 / 32 * [0; 0; 2]; % Komponenten des Kraftvektors anpassen
    neumannObject.meshObject.edof = bounEdof.SY1;

    %% solver
    dofObject = runNewton(setupObject, dofObject);
    centralDisplacement = max(abs(plateObject.qN1(:, 1)));
    disp(['Central displacement: ', num2str(centralDisplacement)]);
    resultArray(ii) = centralDisplacement;

end
