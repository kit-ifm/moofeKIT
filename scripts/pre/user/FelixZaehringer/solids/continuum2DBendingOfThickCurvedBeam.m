% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'huangEtAl2020'}; % {'disp', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
%% Script
for ii = 1:size(elementsToTest, 2)
    innerRadius = 10;
    outerRadius = 15;
    numberOfElements = 4;

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
    setupObject.plotObject.stress.component = -1;
    setupObject.plotObject.view = 2;
    setupObject.usePreconditioning = false;
    setupObject.newton.tolerance = 1e-4;
    setupObject.integrator = 'Endpoint';

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    solidObject = solidClass(dofObject);

    numberOfGausspoints = @(order) (order + 1)^2;

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
        case 'dispPetrovGalerkinCen'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'displacement';
            solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinCenEtAl2015';
        case 'dispPetrovGalerkinXie'
            order = 2;
            serendipity = true;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
            solidObject.mixedFEObject.condensation = true;
            solidObject.mixedFEObject.typeShapeFunctionData = 4;
        case 'easPetrovGalerkin'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinPfefferkornBetsch2021';
            solidObject.elementDisplacementType = 'eas';
            solidObject.mixedFEObject.condensation = false;
            solidObject.mixedFEObject.typeShapeFunctionData = 4;
        case 'huangEtAl2020'
            order = 1;
            serendipity = false;
            solidObject.materialObject.name = 'HookeESZ';
            solidObject.elementDisplacementType = 'incompatibleModes';
            solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinHuangEtAl2020';
            solidObject.mixedFEObject.typeShapeFunctionData = 10;
            numberOfGausspoints = @(order) (order + 2)^2;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    [solidObject.meshObject.nodes, solidObject.meshObject.edof, bounEdof] = meshCurvedBeam(innerRadius, outerRadius, numberOfElements, order, serendipity);

    %plateObject.materialObject.name = 'Hooke';
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.E = 1e3;
    solidObject.materialObject.nu = 0.0;
    solidObject.dimension = 2;
    solidObject.shapeFunctionObject.order = order;
    solidObject.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints(order);

    % Dirichlet boundary
    dirichletBoundary = dirichletClass(dofObject);
    dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
    dirichletBoundary.nodalDof = 1;
    dirichletBoundary.timeFunction = @(t, XYZ) XYZ;
    dirichletBoundary.masterObject = solidObject;

    dirichletBoundary2 = dirichletClass(dofObject);
    dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:, 2) == 0);
    dirichletBoundary2.nodalDof = 2;
    dirichletBoundary2.masterObject = solidObject;

    % Neumann boundary
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = solidObject;
    neumannObject.loadVectorFunction = @(XYZ) [0; 6/5]; % Komponenten des Kraftvektors anpassen
    neumannObject.meshObject.edof = bounEdof.SX1;

    %% solver
    dofObject = runNewton(setupObject, dofObject);
    evalPointX = 0;
    evalPointY = innerRadius;
    nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8);
    endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
    disp(['End displacement: ', num2str(endDisplacement)]);

end
