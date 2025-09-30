allElements = getElementList();
runall = true;

if runall == true
    elements = 1:10;
else
    elements = 9;
end

yDisplacementUpperRightNode = zeros(numel(elements), 1);
for elem = elements
    fprintf('---------------------------------------\n')
    fprintf(['- Current Element: ', allElements(elem).NAME, '\n'])
    fprintf('---------------------------------------\n')

    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'cooksMembrane';
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = 1;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = true;
    setupObject.plotObject.view = 2;
    setupObject.newton.tolerance = 1e-6;
    setupObject.newton.maximumSteps = int32(10);

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    solidObject = solidClass(dofObject);
    solidObject.dimension = 2;
    solidObject.shapeFunctionObject.order = allElements(elem).order;
    solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;
    solidObject.elementDisplacementType = allElements(elem).designator;
    solidObject.mixedFEObject.typeShapeFunctionData = allElements(elem).ansatzIntDOF; %order of pressure DOFs
    solidObject.materialObject.name = 'HookeEVZ';
    solidObject.mixedFEObject.condensation = allElements(elem).condensation;
    solidObject.mixedFEObject.continuousShapeFunctions = allElements(elem).continuousShapeFunctions;
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.E = 1e2;
    solidObject.materialObject.nu = 0.499;
    numberOfElements = 10 / solidObject.shapeFunctionObject.order;
    [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshCooksMembrane(numberOfElements, numberOfElements, solidObject.shapeFunctionObject.order, false);

    % dirchlet boundary conditions
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.masterObject = solidObject;
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = 1:2;
    dirichletObject.timeFunction = str2func('@(t) 0');

    % neumann boundary conditions
    neumannObject = neumannClass(dofObject);
    neumannObject.loadGeometry = 'line';
    neumannObject.masterObject = solidObject;
    neumannObject.loadVector = [0; 5];
    neumannObject.timeFunction = @(t) t;
    neumannObject.meshObject.edof = edofBoundary;

    %% solver
    dofObject = runNewton(setupObject, dofObject);

    %% Postprocessing
    yDisplacementUpperRightNode(elem) = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
    fprintf('\nDisplacement: %4.3f\n\n', yDisplacementUpperRightNode(elem))
end

%#####################################################################################################################################

%% SUBFUNCTIONS
%#####################################################################################################################################
function allElements = getElementList()
% list with all elements
allElements(1:8) = struct('NAME', [], 'designator', [], 'order', 1, 'ansatzIntDOF', 0, 'condensation', true, 'continuousShapeFunctions', false);

%Q1
indx = 1;
allElements(indx).NAME = 'Q1';
allElements(indx).designator = 'displacement';
allElements(indx).order = 1;

%Q2
indx = indx + 1;
allElements(indx).NAME = 'Q2';
allElements(indx).designator = 'displacement';
allElements(indx).order = 2;

%Q1Theta0P0 - full system
indx = indx + 1;
allElements(indx).NAME = 'Q1Theta0P0full';
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 0;
allElements(indx).condensation = false;

%Q1Theta0P0 - static condensation
indx = indx + 1;
allElements(indx).NAME = 'Q1Theta0P0condesation';
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 0;
allElements(indx).condensation = true;

%Q1Theta1P1 - full system, discontinuous
indx = indx + 1;
allElements(indx).NAME = 'Q1Theta1P1discont';
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 1;
allElements(indx).condensation = false;

%Q1Theta1P1 - full system, continuous
indx = indx + 1;
allElements(indx) = allElements(indx-1);
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 1;
allElements(indx).condensation = false;
allElements(indx).continuousShapeFunctions = true;

%Q2Theta1P1 - full system, discontinuous
indx = indx + 1;
allElements(indx).NAME = 'Q2Theta1P1discont';
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 2;
allElements(indx).ansatzIntDOF = 1;
allElements(indx).condensation = false;
allElements(indx).continuousShapeFunctions = false;

%Q2Theta1P1 - full system, continuous
indx = indx + 1;
allElements(indx).NAME = 'Q2Theta1P1cont';
allElements(indx).designator = 'incompressibleSimoTaylorPistor';
allElements(indx).order = 2;
allElements(indx).ansatzIntDOF = 1;
allElements(indx).condensation = false;
allElements(indx).continuousShapeFunctions = true;

%EAS
indx = indx + 1;
allElements(indx).NAME = 'EAS';
allElements(indx).designator = 'eas';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 4;
allElements(indx).condensation = false;
allElements(indx).continuousShapeFunctions = false;

%PianSumihara
indx = indx + 1;
allElements(indx).NAME = 'PianSumihara';
allElements(indx).designator = 'pianSumihara';
allElements(indx).order = 1;
allElements(indx).ansatzIntDOF = 5;
allElements(indx).condensation = true;
allElements(indx).continuousShapeFunctions = false;
end
