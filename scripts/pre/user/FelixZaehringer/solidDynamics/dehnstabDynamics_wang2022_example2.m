%% Parameters
lengthX = 20;
lengthY = 2;
s=1; % 11

% class = 'solid'; % solid, solidVelocity
class = 'solidVelocity'; % solid, solidVelocity

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'dehnstabDynamics';
setupObject.saveObject.saveData = true;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.view = 2;
setupObject.plotObject.postPlotType = 'zero';
setupObject.plotObject.stress.component = 11;
setupObject.newton.tolerance = 1e-2;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'ExplicitEuler';
% setupObject.integrator = 'Newmark';
% setupObject.newmarkGamma         = 0.5;
% setupObject.newmarkBeta          = 0.25;
switch setupObject.integrator
    case {'Endpoint', 'Newmark', 'Midpoint'}
        setupObject.totalTimeSteps = 400;
    case {'ExplicitEuler'}
        setupObject.totalTimeSteps = 10000;
        %courant: dt<=4.5e-4 (12000 steps)
    otherwise
        error('not implemented')
end

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
if strcmpi(class, 'solid')
    solidObject = solidClass(dofObject);
elseif strcmpi(class, 'solidVelocity')
    solidObject = solidVelocityClass(dofObject);
    setupObject.disableTimeIntegration = true;
end
solidObject.dimension = 2;
% solidObject.elementDisplacementType = 'mixed';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
% solidObject.elementNameAdditionalSpecification = 'NewConstraint';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinNewConstraint';
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = false;
% solidObject.mixedFEObject.typeShapeFunction = 1;
% solidObject.mixedFEObject.typeShapeFunctionData = 16;
solidObject.mixedFEObject.condensation = false;
order = 2;
serendipity = true;
solidObject.shapeFunctionObject.order = order;
solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
solidObject.materialObject.name = 'HookeESZ';
if strcmpi(class, 'solid')
    solidObject.materialObject.rho = 10;
elseif strcmpi(class, 'solidVelocity')
    solidObject.materialObject.rho = 100;
    solidObject.materialObject.rhoNew = 100;
end

solidObject.materialObject.E = 1e6;
solidObject.materialObject.nu = 0;

meshType = 'simple';

if strcmpi(meshType, 'simple')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 10, 1, order, serendipity);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
elseif strcmpi(meshType, 'distorted')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 4, 1, order, serendipity);
    nodesToDistortTop = find(nodes(:, 1) ~= -lengthX/2 & nodes(:, 1) ~= lengthX/2 & nodes(:,2) == lengthY/2);
    nodesToDistortBottom = find(nodes(:, 1) ~= -lengthX/2 & nodes(:, 1) ~= lengthX/2 & nodes(:,2) == -lengthY/2);
    nodes(nodesToDistortTop(2:4:end), 1) = nodes(nodesToDistortTop(2:4:end), 1) + s;
    nodes(nodesToDistortTop(4:4:end), 1) = nodes(nodesToDistortTop(4:4:end), 1) - s;
    nodes(nodesToDistortBottom(2:4:end), 1) = nodes(nodesToDistortBottom(2:4:end), 1) - s;
    nodes(nodesToDistortBottom(4:4:end), 1) = nodes(nodesToDistortBottom(4:4:end), 1) + s;
    nodes = adaptMeshToCornerNodes(nodes, solidObject.meshObject.edof, 2);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
end

if strcmpi(class, 'solidVelocity')
    solidObject.meshObject.nodes = [solidObject.meshObject.nodes, [zeros(size(solidObject.meshObject.nodes, 1), 1), zeros(size(solidObject.meshObject.nodes, 1), 1)]];
end

% dirchlet boundary conditions
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.masterObject = solidObject;
dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject1.timeFunction = @(t, XYZ) XYZ;
dirichletObject1.nodalDof = 1;

dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject.timeFunction = @(t, XYZ) XYZ;
dirichletObject.nodalDof = 2;

if strcmpi(class, 'solidVelocity')
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.masterObject = solidObject;
    dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
    dirichletObject1.nodalDof = 3;

    dirichletObject = dirichletClass(dofObject);
    dirichletObject.masterObject = solidObject;
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == lengthX & solidObject.meshObject.nodes(:, 2) == lengthY/2);
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = 4;
end

% neumann boundary
omega = 80;
A = 1;
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [solidObject.materialObject.E * A; 0];
neumannObject.timeFunction = @(t) sin(omega*t);
neumannObject.meshObject.edof = edofBoundary.SX2;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

% FEM solution
nodeToTrack = find(solidObject.meshObject.nodes(:, 1) == lengthX/2 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
xDisplacementNodeToTrack = solidObject.qN1(nodeToTrack, 1) - solidObject.qR(nodeToTrack, 1);

% exact solution
xEval = lengthX/2;
tEval=0:0.001:setupObject.totalTime;
syms x nSym;
c = sqrt(solidObject.materialObject.E/solidObject.materialObject.rho);
p = A*c/(omega*cos(omega*lengthX/c))*sin(omega*x/c);
q_n = int(p*sin((2*nSym-1)/(2*lengthX)*pi*x), x, [0, lengthX]);
xAna = sin(omega*tEval)*eval(subs(p, x, xEval));

for n=1:1000
    A_n = 4*omega*eval(subs(q_n, nSym, n))/((2*n-1)*pi*c);
    xAna = xAna - A_n*sin((2*n-1)/(2*lengthX)*pi*c*tEval)*sin((2*n-1)/(2*lengthX)*pi*xEval);
end

figure;
plot(tEval, xAna);
xlim([0, setupObject.totalTime]);
ylim([-8, 8]);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject, setupObject);
if strcmpi(class, 'solid')
    kineticEnergy = getKineticEnergy(dofObject.postDataObject, setupObject);
elseif strcmpi(class, 'solidVelocity')
    kineticEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'kineticEnergy');
end
strainEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'strainEnergy');
externalEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'externalEnergy');
figure;
plot(timeVector, kineticEnergy+strainEnergy);

%% postprocessing
[setupObjectCell, dofObjectCell] = loadObjectsFromMatFile('dehnstabDynamics', 'all');
nodeToTrack = find(solidObject.meshObject.nodes(:, 1) == lengthX/2 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
xDisplacementNodeToTrack = zeros(size(dofObjectCell, 1), 1);
timeVector = zeros(size(dofObjectCell, 1), 1);
for i=1:size(dofObjectCell, 1)
    setupObjectOfTimestep = setupObjectCell{i};
    dofObjectOfTimestep = dofObjectCell{i};
    solidObjectOfTimestep = dofObjectOfTimestep;
    xDisplacementNodeToTrack(i) = solidObjectOfTimestep(nodeToTrack, 1) - solidObject.qR(nodeToTrack, 1);
    timeVector(i) = setupObjectOfTimestep;
end
figure;
plot([0; timeVector], [0; xDisplacementNodeToTrack]);
xlim([0, setupObject.totalTime]);
ylim([-8, 8]);

