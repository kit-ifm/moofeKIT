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
% Moritz Hille 07.08.2023

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeElasticity';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 800;
setupObject.totalTime = 20;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-8;
setupObject.newton.maximumSteps = 500;

% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
solidObject = solidClass(dofObject);% initialize solidObject
solidObject.meshObject.nodes = abaqusMeshData.qR;
solidObject.meshObject.edof = abaqusMeshData.edof;

% material
% solidViscoObject.materialObject.name = 'Hooke';
solidObject.materialObject.name = 'NeoHooke';
% solidObject.materialObject.name = 'MooneyRivlin';
% solidObject.materialObject.name = 'MooneyRivlinWedge';
%     solidObject.materialObject.name = 'SaintVenant';
%     solidObject.materialObject.name = 'SaintVenantNumericalTangent';
%     solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
solidObject.elementDisplacementType = 'displacement';
%     solidObject.elementDisplacementType = 'eas';
%     solidObject.mixedFEObject.condensation = true;
solidObject.materialObject.rho = 1;
E       = 2.26115e+03;
nu      = 3.95469e-01;
% G (second lame parameter)
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidObject.materialObject.mu = E/(2*(1+nu));

solidObject.materialObject.rho = 10;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;
% solidViscoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidObject.numericalTangentObject.computeNumericalTangent = true;
solidObject.numericalTangentObject.type = 'complex';
% solidObject.numericalTangentObject.showDifferences = true;

bcTimeEnd = 5;
FA = [256;512;768];
FA = FA/10;
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.masterObject = solidObject;
neumannObject1.loadGeometry = 'area';
neumannObject1.loadVector = 1/9*FA;
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.masterObject = solidObject;
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
internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');

% % viscous "energy" -> sum of dissipated energy
% viscousEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'viscousEnergy');
% Dt = setupObject.totalTime/setupObject.totalTimeSteps;
% viscousEnergy = Dt*viscousEnergy;
% for t = 2: setupObject.totalTimeSteps +1
% viscousEnergy(t) = viscousEnergy(t) + viscousEnergy(t-1);
% end

externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = internalEnergy + kineticEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+3:end) - totalEnergy(tStartDiff+2:end-1);
figure;
plot(timeVector, internalEnergy, timeVector, kineticEnergy, timeVector, totalEnergy);
title('internal, kinetic, total elastic, viscous and total energy')
legend('internal energy', 'kinetic energy', 'total elastic energy')
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