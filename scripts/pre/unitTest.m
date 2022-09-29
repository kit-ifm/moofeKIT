%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START SCRIPT FOR UNIT TESTING %%
%%         DO NOT CHANGE!        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% script                    unitTest
% geometry                  L-shaped body
% physics                   continuum-electro-thermo-mechanics
% creator(s)                Felix Zähringer
% date                      2021/12/16
% see also                  LShape, cooksMembrane
 
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'unitTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 6;
setupObject.totalTime = 3;
setupObject.newton.tolerance = 1e-4;
setupObject.plotObject.flag = false;
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityCoarse.inp');
solidElectroThermoObject = solidElectroThermoClass(dofObject);
solidElectroThermoObject.meshObject.nodes = abaqusMeshData.qR;
nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);
nodelistoben = abaqusMeshData.subsets(3).edof(:);
nodelistunten = abaqusMeshData.subsets(4).edof(:);
sizenodesoben = size(nodelistoben(:,1),1);
for l = 1:sizenodesoben
    initialThermalField(nodelistoben(l))=300;
    initialThermalField(nodelistunten(l))=250;
end

initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.meshObject.edof = abaqusMeshData.edof;
solidElectroThermoObject.materialObject.name = 'MooneyRivlin';
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
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = 1;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = true;
% solidElectroThermoObject.numericalTangentObject.type = 'complex';

bcTimeEnd = 2;
FA = [256;512;768];

neumannObject1 = neumannClass(dofObject);
neumannObject1.typeOfLoad = 'deadLoad';
neumannObject1.masterObject = solidElectroThermoObject;
neumannObject1.forceVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;

neumannObject2 = neumannClass(dofObject);
neumannObject2.typeOfLoad = 'deadLoad';
neumannObject2.masterObject = solidElectroThermoObject;
neumannObject2.forceVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

% electrical dirichlet boundary conditions
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 0);
dirichletObject1.nodalDof = 4;
dirichletObject1.masterObject = solidElectroThermoObject;
dirichletObject1.timeFunction = @(t,Z) (6e6).*(t > bcTimeEnd) + (6e6*sin(pi/2*t/bcTimeEnd)).*(t >= 0).*(t <= bcTimeEnd);

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 3);
dirichletObject1.nodalDof = 4;
dirichletObject1.masterObject = solidElectroThermoObject;

%% solver
warning off;
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
internalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
totalEnergy = kineticEnergy + internalEnergy;
lastTotalEnergyValue = totalEnergy(end);
