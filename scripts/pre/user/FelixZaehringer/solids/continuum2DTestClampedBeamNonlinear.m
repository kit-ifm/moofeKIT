% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'pianSumihara'};

%% Script
resultArray = zeros(size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)

    disp('=======================================');
    disp(['Current element: ', elementsToTest{ii}]);
    disp('=======================================');

    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'continuum2DTestClampedBeamNonlinear';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = true;
    setupObject.plotObject.view = 2;
    setupObject.plotObject.postPlotType = 'stress';
    setupObject.plotObject.stress.component = 11;
    setupObject.usePreconditioning = false;
    setupObject.integrator = 'Endpoint';

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    t = 0.05;
    lengthX = 10;
    lengthY = t;
    solidObject = solidClass(dofObject);

    solidObject.numericalTangentObject.computeNumericalTangent = false;
%     solidObject.numericalTangentObject.type = 'centralDifferences';
    solidObject.numericalTangentObject.type = 'complex';
    solidObject.numericalTangentObject.showDifferences = false;

    switch elementsToTest{ii}
        case 'dispQ9'
            order = 2;
            serendipity = false;
            solidObject.elementDisplacementType = 'displacementSC';
            solidObject.materialObject.name = 'SaintVenantEVZ';
        case 'pianSumihara'
            order = 1;
            serendipity = false;
            solidObject.elementDisplacementType = 'pianSumihara';
            solidObject.materialObject.name = 'SaintVenantEVZ';
            solidObject.mixedFEObject.condensation = false;
            solidObject.mixedFEObject.typeShapeFunctionData = 5;
        otherwise
            warning(['Element ', elementsToTest{ii}, ' not implemented!'])
    end

    [solidObject.meshObject.nodes, solidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 10, 1, order, serendipity);
    solidObject.meshObject.nodes(:, 1) = solidObject.meshObject.nodes(:, 1) + lengthX / 2;
    solidObject.meshObject.nodes(:, 2) = solidObject.meshObject.nodes(:, 2) + lengthY / 2;

    solidObject.materialObject.rho = 0;
    E = 1e3;
    nu = 0;
    solidObject.materialObject.lambda = nu/(1-2*nu)*1/(1+nu)*E;
    solidObject.materialObject.mu = 1/2*1/(1+nu)*E;
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
    F = E*t^3/lengthX^3;
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = solidObject;
    neumannObject.loadVector = -[0; F/t]; % Komponenten des Kraftvektors anpassen
    neumannObject.meshObject.edof = bounEdof.SX2;

    %% solver
    dofObject = runNewton(setupObject, dofObject);
    nodeToMeasure = find(solidObject.meshObject.nodes(:, 1) == lengthX & solidObject.meshObject.nodes(:, 2) == 0);
    endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
    disp(['End displacement: ', num2str(endDisplacement)]);

end