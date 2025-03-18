% LSHAPEELECTROTHERMOMECHANICS Script for preprocessing a dynamic coupled 
% electro-thermo-mechanical simulation.
% 
% FORMULATION
% Mixed finite element formulation for nonlinear electro-thermo-mechanical 
% processes is pursued. This is accomplished under the assumption of a 
% hyperelastic, isotropic material and based on a polyconvex inspired 
% internal energy function. The discretization in time is pursued with an 
% energy and momentum conserving scheme. 
% 
% 
% REFERENCE
% https://doi.org/10.1016/j.cma.2021.114298
% https://doi.org/10.1007/BF00913408
% 
% SEE ALSO
% LShape, 
% cooksMembrane
% 
% CREATOR(S) 
% Marlon Franke
 
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeElectroThermo';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 50;
setupObject.totalTime = 40;
setupObject.newton.tolerance = 1e-4;
setupObject.newton.enforceIteration = true;
setupObject.plotObject.flag = false;
% setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.postPlotType = 'phi';
setupObject.plotObject.postPlotType = 'stress';
% setupObject.plotObject.postPlotType = 'D1';
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
% abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
solidElectroThermoObject = solidElectroThermoClass(dofObject);
solidElectroThermoObject.meshObject.nodes = abaqusMeshData.qR;
nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);
% nodelistoben = abaqusMeshData.SUBSETS(3).EDOF(:);
% nodelistunten = abaqusMeshData.SUBSETS(4).EDOF(:);
% sizenodesoben = size(nodelistoben(:,1),1);
% for l = 1:sizenodesoben
%     initialThermalField(nodelistoben(l))=300;
%     initialThermalField(nodelistunten(l))=250;
% end

initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.meshObject.edof = abaqusMeshData.edof;
solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
% solidElectroThermoObject.materialObject.name = 'MooneyRivlin';
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.rho = 100;                                          % mass-density
solidElectroThermoObject.materialObject.a = 831.25;                                         % a
solidElectroThermoObject.materialObject.b = 166.25;                                         % b
solidElectroThermoObject.materialObject.c1 = 10000;                                         % c
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
solidElectroThermoObject.materialObject.thetaR = 293.15;                                  % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 10;                                          % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
% solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.condensation = false;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = 1;
% solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
% solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = 0;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = false;

bcTimeEnd = 5;
FA = [256;512;768];

neumannObject1 = neumannClass(dofObject);
neumannObject1.loadPhysics = 'mechanical';
neumannObject1.loadType = 'deadLoad';
neumannObject1.loadGeometry = 'area';
neumannObject1.masterObject = solidElectroThermoObject;
neumannObject1.loadVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;

neumannObject2 = neumannClass(dofObject);
neumannObject2.loadType = 'deadLoad';
neumannObject2.loadGeometry = 'area';
neumannObject2.masterObject = solidElectroThermoObject;
neumannObject2.loadVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

% %% Test to be removed
% % mechanical dirichlet boundary condition, DirichletFixed
% dirichletObject4 = dirichletClass(dofObject);
% dirichletObject4.nodeList = unique(abaqusMeshData.SUBSETS(4).EDOF(:));
% dirichletObject4.nodalDof = 1:3;
% dirichletObject4.masterObject = solidElectroThermoObject;

% electrical dirichlet boundary conditions
dirichletObject0 = dirichletClass(dofObject);
dirichletObject0.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 0);
dirichletObject0.nodalDof = 4;
dirichletObject0.masterObject = solidElectroThermoObject;
dirichletObject0.timeFunction = str2func('@(t,Z) (6e6).*(t > 5) + (6e6*sin(pi/2*t/5)).*(t >= 0).*(t <= 5)');

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 1.5);
dirichletObject1.nodalDof = 4;
dirichletObject1.masterObject = solidElectroThermoObject;

%% solver
warning off;
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
internalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
DE = getEnergy(dofObject.postDataObject,dofObject,setupObject,'DE');
TS = getEnergy(dofObject.postDataObject,dofObject,setupObject,'TS');
externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
figure; 
%plot(timeVector,kineticEnergy + internalEnergy + DE + TS);
totalEnergy = kineticEnergy + internalEnergy;
plot(timeVector,totalEnergy);

[linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
% figure;
% plot(timeVector,linearMomentum);
figure;
plot(timeVector,totalLinearMomentum);
% figure;
% plot(timeVector,angularMomentum);
figure;
plot(timeVector,totalAngularMomentum);

% matlab2tikz(['diagram','.tikz'],'width','\fwidth','height','\fheight')
