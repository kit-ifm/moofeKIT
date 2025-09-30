% LSHAPE Script for preprocessing a dynamic mechanical simulation.
%
% FORMULATION
% Different formulations like standard displacement-based and mixed en-
% hanced assumed strain (eas) and different material models can be chosen.
%
% REFERENCE
% https://doi.org/10.1007/BF00913408
%
% SEE ALSO
% cooksMembrane,
% LShapeElectroThermo
%
% CREATOR(S)
% Moritz Hille 15.05.2024

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeThermoMechanicsElasticity';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 80;
setupObject.totalTime = 20;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-8;
setupObject.newton.maximumSteps = 500;

setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
solidThermoObject = solidClass(dofObject);% initialize solidObject
% solidThermoObject = solidThermoClass(dofObject);
solidThermoObject.meshObject.nodes = abaqusMeshData.qR;
solidThermoObject.meshObject.edof = abaqusMeshData.edof;

% nnodes=size(solidThermoObject.meshObject.nodes(:,1),1);
% initialThermalField = 293.15*ones(nnodes,1);
% solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, initialThermalField];

% material
% solidViscoObject.materialObject.name = 'Hooke';
% solidThermoObject.materialObject.name = 'NeoHookeGENERIC';
solidThermoObject.materialObject.name = 'MooneyRivlin';
% solidThermoObject.materialObject.name = 'MooneyRivlinWedge';
%     solidThermoObject.materialObject.name = 'SaintVenant';
%     solidThermoObject.materialObject.name = 'SaintVenantNumericalTangent';
%     solidThermoObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
% solidThermoObject.elementDisplacementType = 'displacement';
solidThermoObject.elementDisplacementType = 'displacementSC';
%     solidThermoObject.elementDisplacementType = 'eas';
%     solidThermoObject.mixedFEObject.condensation = true;
solidThermoObject.materialObject.rho = 100;                                                                             % mass density
solidThermoObject.materialObject.a = 1;                                                                                 % material parameter a
solidThermoObject.materialObject.b = 1;                                                                                 % material parameter b
solidThermoObject.materialObject.c = 1;                                                                                 % material parameter c
solidThermoObject.materialObject.d = 2*(solidThermoObject.materialObject.a + 2*solidThermoObject.materialObject.b);     % material parameter d
solidThermoObject.materialObject.kappa = 1;                                                                             % heat capacity

solidThermoObject.materialObject.beta = 0;                                                                  % coupling parameter
% solidThermoObject.materialObject.beta = 2.233*10^(-4);                                                                  % coupling parameter

solidThermoObject.materialObject.thetaR = 293.15;                                                                       % reference Temperature
solidThermoObject.materialObject.k0 = 0.1;                                                                              % thermal conductivity
E       = 2.26115e+03;
nu      = 3.95469e-01;
% G (second lame parameter)
solidThermoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidThermoObject.materialObject.mu = 0;%E/(2*(1+nu));

solidThermoObject.dimension = 3;
solidThermoObject.shapeFunctionObject.order = 1;
solidThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
% solidViscoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidThermoObject.numericalTangentObject.type = 'complex';
% solidThermoObject.numericalTangentObject.showDifferences = true;

bcTimeEnd = 5;
FA = [256;512;768];
FA = FA/10;
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.masterObject = solidThermoObject;
neumannObject1.loadGeometry = 'area';
neumannObject1.loadVector = 1/9*FA;
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.masterObject = solidThermoObject;
neumannObject2.loadGeometry = 'area';
neumannObject2.loadVector = -1/9*FA;
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

%% Test to be removed
% mechanical dirichlet boundary condition, DirichletFixed
% dirichletObject = dirichletClass(dofObject);
% dirichletObject.nodeList = unique(abaqusMeshData.edof);
% dirichletObject.nodalDof = 4;
% dirichletObject.masterObject = solidThermoObject;

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
[linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');

externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = internalEnergy + kineticEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+3:end) - totalEnergy(tStartDiff+2:end-1);
figure;
plot(timeVector, internalEnergy, timeVector, kineticEnergy, timeVector, totalEnergy);
title('internal, kinetic and total energy')
legend('internal energy', 'kinetic energy', 'total energy')
figure;
plot(timeVector(tStartDiff+3:end), totalEnergyDiff);
legend('total energy diff')
title('total energy diff')
matlab2tikz(['diagram','.tikz'],'width','\figW','height','\figH')
figure;
plot(timeVector,linearMomentum);
title('linear momentum')
figure;
plot(timeVector,totalLinearMomentum);
legend('total linear momentum')
title('total linear momentum')
figure;
plot(timeVector,angularMomentum);
title('angular momentum')
figure;
plot(timeVector,totalAngularMomentum);
legend('total angular momentum')
title('total angular momentum')