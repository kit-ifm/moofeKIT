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
% Marlon Franke
%
% Update Visco
% Moritz Hille 25.01.2023

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeViscoelasticity';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 20;
setupObject.totalTime = 4;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-4;
setupObject.newton.maximumSteps = 500;
setupObject.newtonVisco.tolerance = 1e-4;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
solidViscoObject = solidViscoClass(dofObject);% initialize solidViscoObject
solidViscoObject.linearity = 'nonlinear';                                   % nonlinear viscoelastic material behavior
solidViscoObject.meshObject.nodes = abaqusMeshData.qR;
solidViscoObject.meshObject.edof = abaqusMeshData.edof;

% material
solidViscoObject.materialObject.name = 'MooneyRivlinVisco';
solidViscoObject.elementDisplacementType = 'displacement';
solidViscoObject.materialObject.EVisco = 227.2;
solidViscoObject.materialObject.muVisco = 49.875;                              % G (second lame parameter)

bulk    = 5209;                                                 % K
shear   = 997.5;                                                % G
mu      = shear;                                                % second lame parameter
a      = 5/6*mu;                                               % a
b      = 1/6*mu;                                               % b
c       = 10000;
solidViscoObject.materialObject.a = a;
solidViscoObject.materialObject.b = b;
solidViscoObject.materialObject.c = c;
solidViscoObject.materialObject.d = 2*(solidViscoObject.materialObject.a + 2*solidViscoObject.materialObject.b);

solidViscoObject.materialObject.rho = 10;
solidViscoObject.dimension = 3;
solidViscoObject.shapeFunctionObject.order = 1;
solidViscoObject.shapeFunctionObject.numberOfGausspoints = 8;

bcTimeEnd = 5;
FA = [256;512;768];
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.loadType = 'deadLoad';
neumannObject1.masterObject = solidViscoObject;
neumannObject1.loadVector = 1/9*FA;
neumannObject1.loadGeometry = 'area';
neumannObject1.shapeFunctionObject.order = solidViscoObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidViscoObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.loadType = 'deadLoad';
neumannObject2.masterObject = solidViscoObject;
neumannObject2.loadVector = -1/9*FA;
neumannObject2.loadGeometry = 'area';
neumannObject2.shapeFunctionObject.order = solidViscoObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidViscoObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

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
totalEnergyDiff = totalEnergy(tStartDiff+1:end) - totalEnergy(tStartDiff:end-1);
figure;
plot(timeVector, totalEnergy);
title('total energy')
legend('total energy')
figure;
plot(timeVector(tStartDiff+1:end), totalEnergyDiff);
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