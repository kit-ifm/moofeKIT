%% Parameters
lengthX = 50;
lengthY = 50;
s=11; % 11

class = 'solidVelocity'; % solid, solidVelocity

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'block';
setupObject.saveObject.saveData = false;
setupObject.totalTime = 10;
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
        setupObject.totalTimeSteps = 30;
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
end
solidObject.dimension = 2;
solidObject.elementDisplacementType = 'eas';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
% solidObject.elementNameAdditionalSpecification = 'NewConstraint';
solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinNewConstraint';
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = false;
solidObject.mixedFEObject.typeShapeFunction = 2;
solidObject.mixedFEObject.typeShapeFunctionData = 4;
solidObject.mixedFEObject.condensation = true;
order = 1;
serendipity = false;
solidObject.shapeFunctionObject.order = order;
solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
solidObject.materialObject.name = 'HookeESZ';
if strcmpi(class, 'solid')
    solidObject.materialObject.rho = 10;
elseif strcmpi(class, 'solidVelocity')
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.rhoNew = 10;
end

solidObject.materialObject.E = 1e4;
solidObject.materialObject.nu = 0.3;

meshType = 'distorted';

if strcmpi(meshType, 'simple')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 1, 1, order, serendipity);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
elseif strcmpi(meshType, 'distorted')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 2, 2, order, serendipity);
    nodeToDistort = find(nodes(:,1) == 0 & nodes(:,2) == 0);
    nodes(nodeToDistort, :) = nodes(nodeToDistort, :) + s;
    nodes = adaptMeshToCornerNodes(nodes, solidObject.meshObject.edof, 2);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
elseif strcmpi(meshType, 'undistorted')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 2, 2, order, serendipity);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
elseif strcmpi(meshType, 'superfine')
    [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 20, 20, order, serendipity);
    solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
    solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
end

if strcmpi(class, 'solidVelocity')
    solidObject.meshObject.nodes = [solidObject.meshObject.nodes, zeros(size(solidObject.meshObject.nodes))];
end


% dirchlet boundary conditions
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.masterObject = solidObject;
dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject1.nodalDof = 1;

dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
dirichletObject.timeFunction = @(t, XYZ) XYZ;
dirichletObject.nodalDof = 2;

if strcmpi(class, 'solidVelocity')
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.masterObject = solidObject;
    dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
    dirichletObject1.nodalDof = 3;

    dirichletObject = dirichletClass(dofObject);
    dirichletObject.masterObject = solidObject;
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
    dirichletObject.nodalDof = 4;
end

% neumann boundary conditions
loadTime = 1;
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [-800; 0];
neumannObject.timeFunction = @(t) sin(pi*t/loadTime) .* (t <= loadTime);
neumannObject.meshObject.edof = edofBoundary.SX2;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)

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
%
% [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',2);
% figure;
% plot(timeVector, totalAngularMomentum);
%
% [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',2);
% figure;
% plot(timeVector, totalLinearMomentum);
