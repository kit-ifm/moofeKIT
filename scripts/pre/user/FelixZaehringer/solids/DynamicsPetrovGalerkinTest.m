formulation = 'Livens';
% formulation = 'Standard';
% formulation = 'StandardPG';

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'block';
setupObject.saveObject.saveData = false;
setupObject.totalTime = 80;
setupObject.plotObject.flag = false;
setupObject.plotObject.view = 2;
setupObject.newton.tolerance = 1e-2;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'ExplicitEuler';
% setupObject.integrator = 'Newmark';
% setupObject.newmarkGamma         = 0.5;
% setupObject.newmarkBeta          = 0.25;
switch setupObject.integrator
    case {'Endpoint', 'Newmark', 'Midpoint'}
        setupObject.totalTimeSteps = 100;
    case {'ExplicitEuler'}
        setupObject.totalTimeSteps = 10000;
        %courant: dt<=4.5e-4 (12000 steps)
    otherwise
        error('not implemented')
end

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
if strcmp(formulation, 'Livens')
    solidObject = solidLivensClass(dofObject);
    solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
else
    solidObject = solidClass(dofObject);
end
if strcmp(formulation, 'StandardPG')
    solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
end
solidObject.dimension = 2;
solidObject.elementDisplacementType = 'displacement';
% solidObject.elementDisplacementType = 'pianSumihara';
% solidObject.elementDisplacementType = 'selectiveReducedIntegration';
% solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
solidObject.shapeFunctionObject.order = 2;
serendipity = true;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order + 1)^2;
solidObject.materialObject.name = 'HookeEVZ';
% solidObject.materialObject.name = 'Hooke2D';
% solidObject.materialObject.name = 'NeoHooke2D';
solidObject.mixedFEObject.condensation = true;
solidObject.materialObject.E = 1e3;
solidObject.materialObject.nu = 0.2;
% [solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(50, 50, 1, 1, solidObject.shapeFunctionObject.order, serendipity);
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofBoundary] = meshPatchTestDistorted2D(50, 50, solidObject.shapeFunctionObject.order, serendipity);
% solidObject.meshObject.nodes(:, 1) = solidObject.meshObject.nodes(:, 1) + 25;
if strcmp(formulation, 'Livens')
    solidObject.materialObject.rhoLivens = 1;
    solidObject.materialObject.rho = 0;
    solidObject.meshObject.nodes = [solidObject.meshObject.nodes, zeros(size(solidObject.meshObject.nodes)), zeros(size(solidObject.meshObject.nodes))];
    solidObject.numericalTangentObject.computeNumericalTangent = false;
    solidObject.numericalTangentObject.showDifferences = false;
else
    solidObject.materialObject.rho = 1;
end

% dirchlet boundary conditions
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.masterObject = solidObject;
dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject1.nodalDof = 1;

dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find((solidObject.meshObject.nodes(:, 1) == 0) & (abs(solidObject.meshObject.nodes(:, 2)) == 0));
dirichletObject.nodalDof = 1:2;

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line';
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [-200; 0];
neumannObject.timeFunction = @(t) sin(pi*t/2) .* (t <= 1);
neumannObject.meshObject.edof = edofBoundary.SX2;

%% solver
dofObject = runNewton(setupObject, dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject, setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject, setupObject);
strainEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'strainEnergy');
externalEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'externalEnergy');
figure;
if strcmp(formulation, 'Livens')
    plot(timeVector, strainEnergy);
else
    plot(timeVector, kineticEnergy+strainEnergy);
end

% [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',2);
% figure;
% plot(timeVector, totalAngularMomentum);
