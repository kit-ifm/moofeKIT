%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'square';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.elementDisplacementType = 'displacement';
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;
numberOfElements = 3;
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshSquare(numberOfElements, solidObject.shapeFunctionObject.order, 0, 1);
solidObject.materialObject.name = 'HookeESZ';
solidObject.materialObject.rho = 0;
solidObject.materialObject.E = 1e2;
solidObject.materialObject.nu = 0.3;
solidObject.mixedFEObject.condensation = true;
solidObject.dimension = 2;

% boundary conditions - dirichlet
boundaryType = 1;
%1...essentielle RBen
%2...lagrange - punktweise
boundaryNodes = find(solidObject.meshObject.nodes(:, 2) == 0);
if boundaryType == 1
    %x direction
    dirichletObjectX = dirichletClass(dofObject);
    dirichletObjectX.masterObject = solidObject;
    dirichletObjectX.timeFunction = @(t, q) 0;
    dirichletObjectX.nodalDof = 1;
    dirichletObjectX.nodeList = 1;
    %y direction
    dirichletObjectY = dirichletClass(dofObject);
    dirichletObjectY.masterObject = solidObject;
    dirichletObjectY.timeFunction = @(t, q) 0;
    dirichletObjectY.nodalDof = 2;
    dirichletObjectY.nodeList = boundaryNodes;
elseif boundaryType == 2
    %x direction
    dirichletObjectX = dirichletLagrangeClass(dofObject);
    dirichletObjectX.masterObject = solidObject;
    dirichletObjectX.timeFunction = @(t) 0;
    dirichletObjectX.nodalDof = 1;
    dirichletObjectX.nodeList = 1;
    %y direction
    dirichletObjectY = dirichletLagrangeClass(dofObject);
    dirichletObjectY.masterObject = solidObject;
    dirichletObjectY.timeFunction = @(t) 0;
    dirichletObjectY.nodalDof = 2;
    dirichletObjectY.nodeList = boundaryNodes;
end

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = solidObject;
neumannObject.loadGeometry = 'line';
neumannObject.loadVector = [0; -solidObject.materialObject.E / 10];
neumannObject.timeFunction = @(t) t;
neumannObject.meshObject.edof = edofBoundary.SY2;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
%lagerkraft
if boundaryType == 1
    lager = dofObject.R(dirichletObjectY.globalNodesDof);
elseif boundaryType == 2
    lager = dirichletObjectY.qN1;
end
updateNodalForces(neumannObject);
nodalForces = neumannObject.nodalForces;
fprintf('\n')
fprintf('Gesamtlast:            F = %4.4f\n', abs(sum(nodalForces(:, 2))))
fprintf('Gesamt Auflagerkraft:  R = %4.4f\n', abs(sum(lager)))
fprintf('Einzelne Lagerkraefte: ')
fprintf('%4.3f ', lager)
fprintf('\n')
