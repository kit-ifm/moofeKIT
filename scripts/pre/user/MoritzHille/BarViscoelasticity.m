% 3D Bar Script for dynamic mechanical simulation.
%
% FORMULATION
% Different formulations like standard displacement-based and mixed
% polyconvex framework and different material models can be chosen.
%
% SEE ALSO
% LShapeViscoelasticity
% threeDimensionalBar (Bea Hummel)
% 
% CREATOR(S)
% Moritz Hille 08.05.2023

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'BarViscoelasticity';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1000;
setupObject.totalTime = 20;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.postPlotType = 'stress';
setupObject.newton.tolerance = 1e-8;
setupObject.newton.maximumSteps = 500;
setupObject.newtonVisco.tolerance = 1e-8;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidViscoObject = solidViscoClass(dofObject);          % initialize solidViscoObject
solidViscoObject.linearity = 'nonlinear';               % nonlinear viscoelastic material behavior

% material
solidViscoObject.materialObject.name = 'NeoHookeVisco';            
solidViscoObject.elementDisplacementType = 'displacement';

% material parameter
solidViscoObject.materialObject.rho = 10;
E = 2.26115e+03;
nu = 3.95469e-01;
solidViscoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidViscoObject.materialObject.mu = E/(2*(1+nu));
solidViscoObject.materialObject.lambdaVisco = 227.2;
solidViscoObject.materialObject.muVisco = 49.875;  

% mesh
nel = 2;
[solidViscoObject.meshObject.nodes,solidViscoObject.meshObject.edof,edofNeumann] = trilinearBar(4*nel,nel,nel,8,2,2);
solidViscoObject.dimension = 3;
solidViscoObject.shapeFunctionObject.order = 1;
solidViscoObject.shapeFunctionObject.numberOfGausspoints = 8;

% dirichlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1:3;
dirichletObject.masterObject = solidViscoObject;
dirichletObject.timeFunction = str2func('@(t) 0');

% neumann boundary conditions
bcTimeEnd = 5;
FA = [10;0;0];
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = solidViscoObject;
neumannObject.loadGeometry = 'area';
neumannObject.loadVector = FA;
neumannObject2.projectionType = 'none';
neumannObject.timeFunction = @(t) sin(pi * t/bcTimeEnd).*(t<=bcTimeEnd);
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidViscoObject,setupObject)

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

%plots
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