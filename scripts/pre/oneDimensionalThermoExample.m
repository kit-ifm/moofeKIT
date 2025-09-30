%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename();
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 200;
setupObject.totalTime = 2;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = true;
setupObject.plotObject.makeMovie = false;
setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.colorBarLimits = [3200,3300];
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidThermoObject = solidThermoClass(dofObject);
numberOfElements = 5;
[solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof,edofNeumann] = meshOneDimensional(1, numberOfElements, 1);
solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, 293.15*ones(size(solidThermoObject.meshObject.nodes,1),1)];

% solidThermoObject.meshObject.nodes(end,2) = solidThermoObject.meshObject.nodes(end,2) + 20;

solidThermoObject.elementDisplacementType = 'displacement';
solidThermoObject.materialObject.name = 'Hooke1D';
% solidThermoObject.materialObject.name = 'Hooke1DHeatEquation';
solidThermoObject.materialObject.rho = 1;
solidThermoObject.materialObject.E = 1;
solidThermoObject.materialObject.mu = 1;
solidThermoObject.materialObject.kappa = 1;
solidThermoObject.materialObject.beta = 2.233*10^(-4);    % coupling parameter
solidThermoObject.materialObject.thetaR = 293.15;         % reference Temperature
solidThermoObject.materialObject.k0 = 0.1;                % thermal conductivity
solidThermoObject.dimension = 1;
solidThermoObject.shapeFunctionObject.order = 1;
solidThermoObject.shapeFunctionObject.numberOfGausspoints = 2;
solidThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidThermoObject.numericalTangentObject.showDifferences = false;
% solidThermoObject.numericalTangentObject.type = 'complex';

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidThermoObject.meshObject.nodes(:,1)==0);
dirichletObject1.nodalDof = 1;
dirichletObject1.masterObject = solidThermoObject;

nodalLoadObject = nodalLoadClass(dofObject);
nodalLoadObject.masterObject = solidThermoObject;
nodalLoadObject.loadVector = 1;
nodalLoadObject.timeFunction = str2func('@(t) t*(t<1)');
nodalLoadObject.nodeList = find(solidThermoObject.meshObject.nodes(1, :) == 1);

% dirichletObject2 = dirichletClass(dofObject);
% dirichletObject2.nodeList = find(solidThermoObject.meshObject.nodes(:,1)==1);
% dirichletObject2.nodalDof = 1;
% dirichletObject2.masterObject = solidThermoObject;
% dirichletObject2.timeFunction = str2func('@(t,X) (X + 0.5*t*(t<1)) ');

nodalLoadObject2 = nodalLoadClass(dofObject);
nodalLoadObject2.masterObject = solidThermoObject;
nodalLoadObject2.loadPhysics = 'thermal';
nodalLoadObject2.nodeList = 1;
nodalLoadObject2.timeFunction = str2func('@(t) t*1*(t<1)');
nodalLoadObject2.loadVector = 10;

% boundaryK2 = dirichletClass(dofObject);
% boundaryK2.masterObject = solidThermoObject;
% boundaryK2.nodeList = 1;
% boundaryK2.nodalDof = 2;
% boundaryK2.timeFunction = str2func('@(t,temp) temp + t*20*(t<1)');

% boundaryK3 = dirichletClass(dofObject);
% boundaryK3.masterObject = solidThermoObject;
% boundaryK3.nodeList = size(solidThermoObject.qN1,1);
% boundaryK3.nodalDof = 2;
% boundaryK3.timeFunction = str2func('@(t,temp) temp + 20*(t<1)');

% boundaryK4 = dirichletClass(dofObject);
% boundaryK4.masterObject = solidThermoObject;
% boundaryK4.nodeList = 1:size(solidThermoObject.meshObject.nodes,1);
% boundaryK4.nodalDof = 2;
% boundaryK4.timeFunction = str2func('@(t,temp) temp');

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
internalEnergyDifference = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergyDifference');
kineticEnergyDifference = getKineticEnergyDifference(dofObject.postDataObject,setupObject);
totalEnergy = kineticEnergy + internalEnergy;
totalEnergyDifference = kineticEnergyDifference + internalEnergyDifference;
figure;
plot(timeVector,kineticEnergy,timeVector,internalEnergy,timeVector,totalEnergy);
legend('kineticEnergy','internalEnergy','totalEnergy')
xlabel('time [s]')
ylabel('energy [J]')
set(gca,'fontsize',16);
if 0
    matlab2tikz([strcat(setupObject.saveObject.fileName,'Energy'),'.tikz'],'width','\figW','height','\figH')
    print(strcat(setupObject.saveObject.fileName,'Energy'),'-dpng')
end
figure;
plot(timeVector((end-1)/2+1:end),totalEnergyDifference((end-1)/2+1:end));
xlabel('time [s]')
ylabel('energy [J]')
set(gca,'fontsize',16);
if 0
    matlab2tikz([strcat(setupObject.saveObject.fileName,'EnergyDifference'),'.tikz'],'width','\figW','height','\figH')
    print(strcat(setupObject.saveObject.fileName,'EnergyDifference'),'-dpng')
end

% plot(solidThermoObject)
% filename = 'execute';
% VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))
