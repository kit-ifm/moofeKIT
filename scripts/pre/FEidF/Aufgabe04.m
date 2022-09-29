% example = 'patchTestRegular';
example = 'patchTestDistorted';
% example = 'cooksMembrane';

plotBubbleModes = true;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = example;
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
% solidObject.elementDisplacementType = 'displacement';
solidObject.elementDisplacementType = 'incompatibleModesWilson';
% solidObject.elementDisplacementType = 'incompatibleModesTaylor';
solidObject.shapeFunctionObject.order = 1;
solidObject.materialObject.name = 'HookeEVZ';
% solidObject.materialObject.name = 'Hooke';
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1e2;
solidObject.materialObject.nu = 0.499;
solidObject.mixedFEObject.condensation = false;
solidObject.dimension = 2;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;

switch example
    case 'patchTestRegular'
        [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshSquare(3, solidObject.shapeFunctionObject.order, 0, 1);

        % dirchlet boundary conditions 1
        dirichletObject1 = dirichletClass(dofObject);
        dirichletObject1.dimension = 2;
        dirichletObject1.masterObject = solidObject;
        dirichletObject1.nodeList = find((solidObject.meshObject.nodes(:, 1) == 0).*(solidObject.meshObject.nodes(:, 2) == 0));
        dirichletObject1.nodalDof = 1;
        dirichletObject1.timeFunction = str2func('@(t,XYZ) 0');

        % dirchlet boundary conditions 2
        dirichletObject2 = dirichletClass(dofObject);
        dirichletObject2.dimension = 2;
        dirichletObject2.masterObject = solidObject;
        dirichletObject2.nodeList = unique(edofBoundary.SY1(:));
        dirichletObject2.nodalDof = 2;
        dirichletObject2.timeFunction = str2func('@(t,XYZ) 0');

        % neumann boundary conditions
        neumannObject = neumannClass(dofObject);
        neumannObject.dimension = 2;
        neumannObject.typeOfLoad = 'deadLoad';
        neumannObject.masterObject = solidObject;
        neumannObject.forceVector = [0; -10];
        neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
        neumannObject.shapeFunctionObject.numberOfGausspoints = neumannObject.shapeFunctionObject.order + 1;
        neumannObject.projectionType = 'none';
        neumannObject.timeFunction = @(t) t;
        neumannObject.meshObject.edof = edofBoundary.SY2;
    case 'patchTestDistorted'
        solidObject.meshObject.nodes = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 6, 6; ...
            0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 0, 3, 0, 3]';
        solidObject.meshObject.edof = [1, 5, 6, 2; ...
            2, 6, 7, 3; ...
            3, 7, 8, 4; ...
            5, 11, 9, 6; ...
            6, 9, 10, 7; ...
            7, 10, 12, 8; ...
            9, 11, 12, 10; ...
            11, 13, 14, 12];
        % dirchlet boundary conditions 1
        dirichletObject1 = dirichletClass(dofObject);
        dirichletObject1.dimension = 2;
        dirichletObject1.masterObject = solidObject;
        dirichletObject1.nodeList = 1;
        dirichletObject1.nodalDof = 2;
        dirichletObject1.timeFunction = str2func('@(t,XYZ) 0');

        % dirchlet boundary conditions 2
        dirichletObject2 = dirichletClass(dofObject);
        dirichletObject2.dimension = 2;
        dirichletObject2.masterObject = solidObject;
        dirichletObject2.nodeList = 1:4;
        dirichletObject2.nodalDof = 1;
        dirichletObject2.timeFunction = str2func('@(t,XYZ) 0');

        %boundary conditions - neumann
        neumannObject = neumannClass(dofObject);
        neumannObject.dimension = 2;
        neumannObject.masterObject = solidObject;
        neumannObject.forceVector = [-30; 0];
        neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
        neumannObject.shapeFunctionObject.numberOfGausspoints = neumannObject.shapeFunctionObject.order + 1;
        neumannObject.projectionType = 'none';
        neumannObject.timeFunction = @(t) t;
        neumannObject.meshObject.edof = [13, 14];
    case 'cooksMembrane'
        numberOfElements = 10;
        [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshCooksMembrane(numberOfElements, numberOfElements, solidObject.shapeFunctionObject.order);

        % dirchlet boundary conditions
        dirichletObject = dirichletClass(dofObject);
        dirichletObject.dimension = 2;
        dirichletObject.masterObject = solidObject;
        dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
        dirichletObject.nodalDof = 1:2;
        dirichletObject.timeFunction = str2func('@(t) 0');

        % neumann boundary conditions
        neumannObject = neumannClass(dofObject);
        neumannObject.dimension = 2;
        neumannObject.typeOfLoad = 'deadLoad';
        neumannObject.masterObject = solidObject;
        neumannObject.forceVector = [0; 5];
        neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
        neumannObject.shapeFunctionObject.numberOfGausspoints = neumannObject.shapeFunctionObject.order + 1;
        neumannObject.projectionType = 'none';
        neumannObject.timeFunction = @(t) t;
        neumannObject.meshObject.edof = edofNeumann;
end

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)