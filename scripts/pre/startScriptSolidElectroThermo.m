% STARTSCRIPTSOLIDELECTROTHERMO Script for preprocessing a minimal dynamic
% coupled electro-thermo-mechanical simulation.
%
% FORMULATION
% Mixed finite element formulation for nonlinear electro-thermo-mechanical
% processes is pursued. This is accomplished under the assumption of a
% hyperelastic, isotropic material and based on a polyconvex inspired
% internal energy function. The discretization in time is pursued with an
% energy and momentum conserving scheme.
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2021.114298
%
% SEE ALSO
% LShapeElectroThermoMechanics
%
% CREATOR(S)
% Marlon Franke

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solidElectroThermo';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 4;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.colorBarLimits = [3200,3300];
setupObject.plotObject.makeMovie = false;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
solidElectroThermoObject = solidElectroThermoClass(dofObject);
[solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1,1,1,1,1,1,2,true);
solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
nnodes = size(solidElectroThermoObject.meshObject.nodes(:, 1), 1);
initialThermalField = 293.15 * ones(nnodes, 1);
initialElectricalPotentialField = zeros(nnodes, 1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
% solidElectroThermoObject.elementDisplacementType = 'mixedSC';
% solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.name = 'MooneyRivlin';
solidElectroThermoObject.materialObject.rho = 0; % mass-density
solidElectroThermoObject.materialObject.a = 831.25; % a
solidElectroThermoObject.materialObject.b = 166.25; % b
solidElectroThermoObject.materialObject.c1 = 10000; % c
solidElectroThermoObject.materialObject.c2 = 5209; % e
solidElectroThermoObject.materialObject.d1 = 2 * (solidElectroThermoObject.materialObject.a + 2 * solidElectroThermoObject.materialObject.b); % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541 * 10^(-12); % epsilon_0
solidElectroThermoObject.materialObject.e1 = 4; % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 100; % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233 * 10^(-4); % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 293.15; % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 10; % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = 1;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = true;
solidElectroThermoObject.numericalTangentObject.showDifferences = true;
solidElectroThermoObject.numericalTangentObject.type = 'complex';

dirichletBoundary1 = dirichletClass(dofObject);
dirichletBoundary1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 3) == 0);
dirichletBoundary1.nodalDof = 3;
dirichletBoundary1.masterObject = solidElectroThermoObject;

dirichletBoundary2 = dirichletClass(dofObject);
dirichletBoundary2.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 1) == 0);
dirichletBoundary2.nodalDof = 1;
dirichletBoundary2.masterObject = solidElectroThermoObject;

dirichletBoundary3 = dirichletClass(dofObject);
dirichletBoundary3.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 2) == 0);
dirichletBoundary3.nodalDof = 2;
dirichletBoundary3.masterObject = solidElectroThermoObject;

dirichletBoundary4 = dirichletClass(dofObject);
dirichletBoundary4.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 3) == 1);
dirichletBoundary4.nodalDof = 3;
dirichletBoundary4.masterObject = solidElectroThermoObject;
dirichletBoundary4.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

dirichletBoundary5 = dirichletClass(dofObject);
dirichletBoundary5.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 2) == 0);
dirichletBoundary5.nodalDof = 4;
dirichletBoundary5.masterObject = solidElectroThermoObject;
dirichletBoundary5.timeFunction = str2func('@(t,Z) 1e6');

dirichletBoundary6 = dirichletClass(dofObject);
dirichletBoundary6.nodeList = find(solidElectroThermoObject.meshObject.nodes(:, 2) == 1);
dirichletBoundary6.nodalDof = 4;
dirichletBoundary6.masterObject = solidElectroThermoObject;

%% solver
warning off
dofObject = runNewton(setupObject, dofObject);

%% postprocessing - energy
% f = figure;
% plot(solidElectroThermoObject, setupObject);
% xlim([0 inf]);
% ylim([0 inf]);
% zlim([0 inf]);
% xticks([]);
% yticks([]);
% zticks([]);
% f.Position = [100 100 320 210];
% zlim([0 1])
% xlabel('x');
% ylabel('y');
% zlabel('z');
% caxis([-1 1]);
% campos([-5,-6,5]);
% % export_fig('test','-eps');
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
internalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = kineticEnergy + internalEnergy;
figure; 
plot(timeVector,totalEnergy);
%plot(solidThermoObject)
% filename = 'execute';
% VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))