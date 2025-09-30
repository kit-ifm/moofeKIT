%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solidThermo';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 20;
setupObject.totalTime = 2;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = true;
setupObject.plotObject.makeMovie = false;
setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.colorBarLimits = [3200,3300];
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidThermoObject = solidThermoClass(dofObject);

numberOfElements = 2;
[solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, edofBoundary] = meshGeneratorCube(1,1,1,2,2,2,1,false);
solidThermoObject.meshObject.nodes = solidThermoObject.meshObject.nodes + 0.5;
solidThermoObject.shapeFunctionObject.numberOfGausspoints = 27;

% [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof] = triquadraticTetrahedralBrick(numberOfElements, numberOfElements, numberOfElements, 1, 1, 1, 1);
% solidThermoObject.meshObject.nodes(:,3) = solidThermoObject.meshObject.nodes(:,3) + 1;
% solidThermoObject.shapeFunctionObject.numberOfGausspoints = 11;

solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, 293.15*ones(size(solidThermoObject.meshObject.nodes,1),1)];
solidThermoObject.elementDisplacementType = 'displacementSC';
solidThermoObject.materialObject.name = 'MooneyRivlin';
solidThermoObject.materialObject.rho = 0;
solidThermoObject.materialObject.a = 1;
solidThermoObject.materialObject.b = 1;
solidThermoObject.materialObject.c = 1;
solidThermoObject.materialObject.d = 2*(solidThermoObject.materialObject.a + 2*solidThermoObject.materialObject.b);
solidThermoObject.materialObject.kappa = 1;
solidThermoObject.materialObject.beta = 2.233*10^(-4);    % coupling parameter
solidThermoObject.materialObject.thetaR = 293.15;         % reference Temperature
solidThermoObject.materialObject.k0 = 0.1;                % thermal conductivity
solidThermoObject.dimension = 3;
solidThermoObject.shapeFunctionObject.order = 2;
solidThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidThermoObject.numericalTangentObject.showDifferences = false;
solidThermoObject.numericalTangentObject.type = 'complex';

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidThermoObject.meshObject.nodes(:,3) == 1);
dirichletObject1.nodalDof = 3;
dirichletObject1.masterObject = solidThermoObject;
dirichletObject1.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = find(solidThermoObject.meshObject.nodes(:,3) == 0);
dirichletObject2.nodalDof = 3;
dirichletObject2.masterObject = solidThermoObject;

dirichletObject3 = dirichletClass(dofObject);
dirichletObject3.nodeList = find(solidThermoObject.meshObject.nodes(:,1) == 0);
dirichletObject3.nodalDof = 1;
dirichletObject3.masterObject = solidThermoObject;

dirichletObject4 = dirichletClass(dofObject);
dirichletObject4.nodeList = find(solidThermoObject.meshObject.nodes(:,2) == 0);
dirichletObject4.nodalDof = 2;
dirichletObject4.masterObject = solidThermoObject;

neumannObject = neumannClass(dofObject);
neumannObject.masterObject = solidThermoObject;
neumannObject.meshObject.edof = edofBoundary.SZ2;
neumannObject.loadVector = 10;
neumannObject.loadPhysics = 'thermal';
neumannObject.loadGeometry = 'area';

% Body force
% bodyLoad = bodyForceClass(dofObject);
% bodyLoad.masterObject = solidThermoObject;
% bodyLoad.meshObject.nodes = solidThermoObject.meshObject.nodes;
% bodyLoad.meshObject.edof = solidThermoObject.meshObject.edof;
% bodyLoad.shapeFunctionObject.order = solidThermoObject.shapeFunctionObject.order;
% bodyLoad.shapeFunctionObject.numberOfGausspoints = solidThermoObject.shapeFunctionObject.numberOfGausspoints;
% bodyLoad.typeOfLoad = 'deadLoad';
% bodyLoad.timeFunction = @(t) t;
% bodyLoad.loadFunction = [0;0;-9.81];
%
% boundaryK2 = dirichletClass(dofObject);
% boundaryK2.masterObject = solidThermoObject;
% boundaryK2.nodeList = 1:size(solidThermoObject.meshObject.nodes,1);
% boundaryK2.nodalDof = 4;
% boundaryK2.timeFunction = str2func('@(t,temp) temp');

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
% timeVector = getTime(dofObject.postDataObject,setupObject);
% kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
% strainEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
% externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
% figure;
% plot(timeVector,kineticEnergy + strainEnergy);
% %plot(solidThermoObject)
% filename = 'execute';
% VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))
