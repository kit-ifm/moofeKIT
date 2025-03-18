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
setupObject.totalTimeSteps = 80;
setupObject.totalTime = 20;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-8;
setupObject.newton.maximumSteps = 500;
setupObject.newtonVisco.tolerance = 1e-8;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
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
% solidViscoObject.materialObject.name = 'Hooke';
solidViscoObject.materialObject.name = 'NeoHookeVisco';
%     solidObject.materialObject.name = 'SaintVenant';
%     solidObject.materialObject.name = 'SaintVenantNumericalTangent';
%     solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
% solidViscoObject.elementDisplacementType = 'displacement'; 
    solidViscoObject.elementDisplacementType = 'displacement';
%     solidObject.elementDisplacementType = 'eas';
%     solidObject.mixedFEObject.condensation = true;
solidViscoObject.materialObject.rho = 1;
E       = 2.26115e+03;
nu      = 3.95469e-01;
% EVisco  = E;
% nuVisco = nu;
solidViscoObject.materialObject.lambdaVisco = 227.2;
solidViscoObject.materialObject.muVisco = 49.875;                              % G (second lame parameter)
% solidViscoObject.materialObject.lambdaVisco = 0;
% solidViscoObject.materialObject.muVisco = 0;
solidViscoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidViscoObject.materialObject.mu = E/(2*(1+nu));

solidViscoObject.materialObject.rho = 10;
solidViscoObject.dimension = 3;
solidViscoObject.shapeFunctionObject.order = 1;
solidViscoObject.shapeFunctionObject.numberOfGausspoints = 8;
% solidViscoObject.shapeFunctionObject.numberOfGausspoints = 27;
% solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.showDifferences = true;

bcTimeEnd = 5;
FA = [256;512;768];
% FA = FA/10;
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.masterObject = solidViscoObject;
neumannObject1.loadGeometry = 'area';
neumannObject1.loadVector = 1/9*FA;
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.masterObject = solidViscoObject;
neumannObject2.loadGeometry = 'area';
neumannObject2.loadVector = -1/9*FA;
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
internalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');

% viscous "energy" -> sum of dissipated energy
viscousEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'viscousEnergy');
Dt = setupObject.totalTime/setupObject.totalTimeSteps;
viscousEnergy = Dt*viscousEnergy;
for t = 2: setupObject.totalTimeSteps +1
viscousEnergy(t) = viscousEnergy(t) + viscousEnergy(t-1);
end

externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
ElasticEnergy = internalEnergy + kineticEnergy;
totalEnergy = ElasticEnergy + viscousEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+1:end) - totalEnergy(tStartDiff:end-1);
figure;
plot(timeVector, internalEnergy, timeVector, kineticEnergy, timeVector, ElasticEnergy, timeVector, viscousEnergy, timeVector, totalEnergy);
title('internal, kinetic, total elastic, viscous and total energy')
legend('internal energy', 'kinetic energy', 'total elastic energy', 'viscous energy', 'total energy')
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