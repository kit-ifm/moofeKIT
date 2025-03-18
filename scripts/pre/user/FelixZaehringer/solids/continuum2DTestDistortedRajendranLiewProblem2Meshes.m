% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'dispPetrovGalerkinCen'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
meshType = 2;

%% Script
for ii = 1:size(elementsToTest, 2)

    disp('=======================================');
    disp(['Current element: ', elementsToTest{ii}]);
    disp('=======================================');

    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'continuum2DTestDistorted';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = true;
    setupObject.plotObject.postPlotType = 'stress';
    setupObject.plotObject.stress.component = 11;
    setupObject.plotObject.view = 2;
    setupObject.usePreconditioning = false;
    setupObject.newton.tolerance = 1e-4;
    setupObject.integrator = 'Endpoint';

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    lengthX = 100;
    lengthY = 10;
    if meshType == 4
        lengthX = 20;
    end
    solidObject = solidClass(dofObject);

    switch elementsToTest{ii}
        case 'disp'
            order = 2;
            serendipity = true;
            solidObject.materialObject.name = 'HookeESZ';
        case 'dispQ4'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
        case 'dispQ9'
            order = 2;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
        case 'eas'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'eas';
            solidObject.mixedFEObject.condensation = true;
            solidObject.mixedFEObject.typeShapeFunctionData = 4;
        case 'dispPetrovGalerkin'
            order = 2;
            serendipity = true;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'displacementPetrovGalerkin';
        case 'dispPetrovGalerkinXie'
            order = 2;
            serendipity = true;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
            solidObject.mixedFEObject.condensation = true;
            solidObject.mixedFEObject.typeShapeFunctionData = 4;
        case 'dispPetrovGalerkinCen'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'displacement';
            solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinCenEtAl2015';
        case 'easPetrovGalerkin'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
            solidObject.elementDisplacementType = 'eas';
            solidObject.mixedFEObject.condensation = false;
            solidObject.mixedFEObject.typeShapeFunctionData = 4;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    [solidObject.meshObject.nodes, solidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 2, 1, order, serendipity);
    if meshType == 1
        % nothing to do here
    elseif meshType == 2
        s = 25;
        nodeToDistort1 = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == -lengthY/2);
        nodeToDistort2 = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
        solidObject.meshObject.nodes(nodeToDistort1, 1) = solidObject.meshObject.nodes(nodeToDistort1, 1) - s;
        solidObject.meshObject.nodes(nodeToDistort2, 1) = solidObject.meshObject.nodes(nodeToDistort2, 1) + s;
        solidObject.meshObject.nodes = adaptMeshToCornerNodes(solidObject.meshObject.nodes, solidObject.meshObject.edof, 2);
    elseif meshType == 4
        s = 2;
        nodeToDistort = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == 0);
        solidObject.meshObject.nodes(nodeToDistort, 1) = solidObject.meshObject.nodes(nodeToDistort, 1) - s;
        if order == 2 && ~serendipity
            solidObject.meshObject.nodes(nodeToDistort - 1, 1) = solidObject.meshObject.nodes(nodeToDistort - 1, 1) - s / 2;
            solidObject.meshObject.nodes(nodeToDistort + 1, 1) = solidObject.meshObject.nodes(nodeToDistort + 1, 1) - s / 2;
        end
    elseif meshType == 5
        s = 7.5;
        nodeToDistort1 = find(solidObject.meshObject.nodes(:, 1) == -lengthX/4 & solidObject.meshObject.nodes(:, 2) == -lengthY/2);
        nodeToDistort2 = find(solidObject.meshObject.nodes(:, 1) == lengthX/4 & solidObject.meshObject.nodes(:, 2) == -lengthY/2);
        solidObject.meshObject.nodes(nodeToDistort1, 1) = solidObject.meshObject.nodes(nodeToDistort1, 1) + s;
        solidObject.meshObject.nodes(nodeToDistort2, 1) = solidObject.meshObject.nodes(nodeToDistort2, 1) + s;
        if order == 2 && ~serendipity
            nodeToDistort3 = find(solidObject.meshObject.nodes(:, 1) == -lengthX/4 & solidObject.meshObject.nodes(:, 2) == 0);
            nodeToDistort4 = find(solidObject.meshObject.nodes(:, 1) == lengthX/4 & solidObject.meshObject.nodes(:, 2) == 0);
            solidObject.meshObject.nodes(nodeToDistort3, 1) = solidObject.meshObject.nodes(nodeToDistort3, 1) + s / 2;
            solidObject.meshObject.nodes(nodeToDistort4, 1) = solidObject.meshObject.nodes(nodeToDistort4, 1) + s / 2;
        end
    end
    solidObject.meshObject.nodes(:, 1) = solidObject.meshObject.nodes(:, 1) + lengthX / 2;
    solidObject.meshObject.nodes(:, 2) = solidObject.meshObject.nodes(:, 2) + lengthY / 2;

    %plateObject.materialObject.name = 'Hooke';
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.E = 1e7;
    solidObject.materialObject.nu = 0.3;
    solidObject.dimension = 2;
    solidObject.shapeFunctionObject.order = order;
    solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

    % Dirichlet boundary
    dirichletBoundary = dirichletClass(dofObject);
    dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == 0);
    dirichletBoundary.nodalDof = 1:2;
    dirichletBoundary.masterObject = solidObject;

    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY);
    dirichletBoundary2.nodalDof = 1;
    dirichletBoundary2.masterObject = solidObject;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = solidObject;
    neumannObject.loadVectorFunction = @(XYZ) [0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
    neumannObject.meshObject.edof = bounEdof.SX2;

    neumannObject2 = neumannClass(dofObject);
    neumannObject2.loadGeometry = 'line';
    neumannObject2.masterObject = solidObject;
    neumannObject2.loadVectorFunction = @(XYZ) -[0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
    neumannObject2.meshObject.edof = bounEdof.SX1;

    neumannObject3 = neumannClass(dofObject);
    neumannObject3.loadGeometry = 'line';
    neumannObject3.masterObject = solidObject;
    neumannObject3.loadVectorFunction = @(XYZ) [240*XYZ(2)/lengthY-120; 0]; % Komponenten des Kraftvektors anpassen
    neumannObject3.meshObject.edof = bounEdof.SX1;

    %% solver
    dofObject = runNewton(setupObject, dofObject);
    pointToMeasureX = lengthX;
    pointToMeasureY = 0;
    nodeToMeasure = find(solidObject.meshObject.nodes(:, 1) == pointToMeasureX & solidObject.meshObject.nodes(:, 2) == pointToMeasureY);
    endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
    exactSolution = (-40/(lengthY*lengthX)*pointToMeasureX^3-36/(lengthY*lengthX)*pointToMeasureX*pointToMeasureY^2+120/lengthY*pointToMeasureX^2+36/lengthX*pointToMeasureX*pointToMeasureY+36/lengthY*pointToMeasureY^2+46*lengthY/lengthX*pointToMeasureX-36*pointToMeasureY)/solidObject.materialObject.E;
    disp(['End displacement: ', num2str(endDisplacement)]);
    disp(['Ratio: ', num2str(endDisplacement/exactSolution)]);
end