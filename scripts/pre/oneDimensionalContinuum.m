%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime = 5;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = 1;
setupObject.newton.tolerance = 1e-4;
% setupObject.integrator = 'Endpoint'; 
setupObject.integrator = 'Midpoint'; 

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
    %für linear-elastisch:
% solidViscoObject = solidClass(dofObject);
    %für linear-viskoelastisch:
solidViscoObject = solidViscoClass(dofObject);
solidViscoObject.linearity = 'linear';
% solidViscoObject.materialObject.name = 'HookeSplit';
solidViscoObject.materialObject.name = 'Hooke';
solidViscoObject.elementDisplacementType = 'displacement';
% E = 2.26115e+03;
 E = 300000;
nu = 3e-01;
solidViscoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));                 % Lame parameter
solidViscoObject.materialObject.mu = E/(2*(1+nu));
% solidViscoObject.materialObject.rho = 10; 
solidViscoObject.materialObject.rho = 1500; 

solidViscoObject.materialObject.eModul0 = 250000;
solidViscoObject.materialObject.eModul1 = 50000;
solidViscoObject.materialObject.eta1 = 10000;

nel = 20;

%Laenge des Stabes = 1
[solidViscoObject.meshObject.nodes,solidViscoObject.meshObject.edof,edofNeumann] = meshOneDimensional(1,nel,1);
solidViscoObject.meshObject.nodes = solidViscoObject.meshObject.nodes + 1/2;
% Warum macht diese Zeile so einen Unterschied? --> wegen dem Dirichlet-Rand, hier nur ein Knoten

solidViscoObject.dimension = 1;
solidViscoObject.shapeFunctionObject.order = 1;
solidViscoObject.shapeFunctionObject.numberOfGausspoints = 2;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1;
dirichletObject.masterObject = solidViscoObject;
dirichletObject.timeFunction = str2func('@(t,X) 0*X');
dirichletObject.dimension = 1;

neumannObject = neumannClass(dofObject);
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidViscoObject;
% neumannObject.shapeFunctionObject.order = solidVscoObject.shapeFunctionObject.order;
% neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidViscoObject.dimension-1);
neumannObject.projectionType = 'none';
neumannObject.forceVector = [10000];
% tol = 1e-12; % numerische Toleranz für if-Bedingung
neumannObject.timeFunction = str2func('@(t)  sin(pi/2 * t).*(t<=1+1e-12)');
neumannObject.meshObject.edof = edofNeumann;
neumannObject.dimension = 1;



% dirichletObject = dirichletClass(dofObject);
% dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==1);
% dirichletObject.nodalDof = 1;
% dirichletObject.masterObject = solidViscoObject;
% dirichletObject.timeFunction = str2func('@(t,X) (X + 0.5*t*(t<1)) ');
% dirichletObject.dimension = 1;

%% solver
dofObject = runNewton(setupObject,dofObject);
plot(solidViscoObject,setupObject)


%% Energien auswerten
kineticEnergy = zeros(setupObject.totalTimeSteps,1);
potentialEnergy = zeros(setupObject.totalTimeSteps,1);
timeStep = zeros(setupObject.totalTimeSteps,1);
dissipatedWork = zeros(setupObject.totalTimeSteps,1);

for j = 1:setupObject.totalTimeSteps
    timeStep(j) = j;
    kineticEnergy(j) = getfield(dofObject.postDataObject.energyJournal(j+1), 'EKin');
    potentialEnergy(j) =  getfield(dofObject.listContinuumObjects{1,1}.ePot(j+1), 'strainEnergy');
    dissipatedWork(j) = getfield(dofObject.listContinuumObjects{1,1}.ePot(j+1), 'dissipatedWork');
end

figure
plot(timeStep, kineticEnergy, timeStep, potentialEnergy, timeStep, kineticEnergy+potentialEnergy, timeStep, dissipatedWork, timeStep, kineticEnergy+potentialEnergy+dissipatedWork)
legend('Ekin', 'Epot', 'Efree', 'Summe Wdiss', 'Eges');

Energiedifferenz = zeros((setupObject.totalTimeSteps-1),1);

f = setupObject.totalTimeSteps/5+1:(setupObject.totalTimeSteps-1);

for i= 1:(setupObject.totalTimeSteps-1)
    Energiedifferenz(i) = (kineticEnergy(i+1)-kineticEnergy(i))+(potentialEnergy(i+1)-potentialEnergy(i))+dissipatedWork(i+1)-dissipatedWork(i); 
end

figure
plot(f,Energiedifferenz(f))
legend('Energiedifferenz');
