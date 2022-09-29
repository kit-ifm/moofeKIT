% ELECTROTHERMOMECHANICSACTUATIONTEST Script for preprocessing a dynamic
% coupled electro-thermo-mechanical simulation.
%
% REFERENCE
% In preparation
%
% SEE ALSO
% LShapeElectroThermoMechanics
%
% CREATOR(S) 
% Felix Zähringer, Marlon Franke

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'electroThermoMechanicsActuationTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime = 1;
setupObject.newton.tolerance = 1e-4;
setupObject.plotObject.flag = true;
setupObject.plotObject.makeMovie = false;
setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.postPlotType = 'phi';
% setupObject.plotObject.postPlotType = 'stress';
% setupObject.plotObject.postPlotType = 'none';
% setupObject.plotObject.postPlotType = 'D1';
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-46,20];
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidElectroThermoObject = solidElectroThermoClass(dofObject);
serendipity = true;
order = 2;
numberOfElementsX = 10;
numberOfElementsY = 10;
numberOfElementsZ = 4;
[solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1, 1, 0.2, numberOfElementsX, numberOfElementsY, numberOfElementsZ, order, serendipity);
solidElectroThermoObject.meshObject.nodes(:,1:2) = solidElectroThermoObject.meshObject.nodes(:,1:2) + 0.5;
solidElectroThermoObject.meshObject.nodes(:,3) = solidElectroThermoObject.meshObject.nodes(:,3) + 0.1;
theta0 = 293.15;
initialThermalField = theta0*ones(size(solidElectroThermoObject.meshObject.nodes,1),1);
initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.rho = 1000;                                         % mass-density
solidElectroThermoObject.materialObject.a = 25000;                                         % a
solidElectroThermoObject.materialObject.b = 50000;                                         % b
solidElectroThermoObject.materialObject.c1 = 500000;                                         % c
solidElectroThermoObject.materialObject.c2 = 5209;                                          % e
solidElectroThermoObject.materialObject.d1 = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12);                             	% epsilon_0
solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 100;                                        % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                             % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 273.15;                                  % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 10;                                          % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.orderShapeFunction = 1;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = false;

% mechanical dirichlet boundary conditions
dirichletObjectMechanical = dirichletClass(dofObject);
dirichletObjectMechanical.masterObject = solidElectroThermoObject;
dirichletObjectMechanical.nodeList = unique([find((round(solidElectroThermoObject.meshObject.nodes(:,3),2) == 0).*(round(solidElectroThermoObject.meshObject.nodes(:,1),2) == 0)),...
                                    find((round(solidElectroThermoObject.meshObject.nodes(:,3),2) == 0).*(round(solidElectroThermoObject.meshObject.nodes(:,1),2) == 1))]);
dirichletObjectMechanical.nodalDof = 1:3;

% electrical dirichlet boundary conditions
dirichletObjectElectrical1 = dirichletClass(dofObject);
dirichletObjectElectrical1.masterObject = solidElectroThermoObject;
dirichletObjectElectrical1.nodeList = find(round(solidElectroThermoObject.meshObject.nodes(:,3),2) == 0);
dirichletObjectElectrical1.nodalDof = 4;

dirichletObjectElectrical2 = dirichletClass(dofObject);
dirichletObjectElectrical2.nodeList = find(round(solidElectroThermoObject.meshObject.nodes(:,3),2) == 0.1);
dirichletObjectElectrical2.nodalDof = 4;
maxPotential = 6e6;
tendLoad = 0.4;
dirichletObjectElectrical2.timeFunction = @(t,Z) (maxPotential*sin(pi/2*t/tendLoad)).*(t >= 0).*(t <= tendLoad);
dirichletObjectElectrical2.masterObject = solidElectroThermoObject;

%% solver
warning off;
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
dofObject = runNewton(setupObject,dofObject);
