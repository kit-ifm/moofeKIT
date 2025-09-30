

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename();
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 200;
setupObject.totalTime = 3600*30;
setupObject.newton.tolerance = 1e-4;
setupObject.plotObject.flag = true;
setupObject.plotObject.makeMovie = false;
setupObject.plotObject.postPlotType = 'temp';
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.colorBarLimits = [280,330];
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidThermoObject = solidThermoClass(dofObject);
nelX = 2;
nelY = 2;
nelZ = 1;
E = 100e+5;
ny = 0.3;
alpha = 1e-3;

[solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nelX ,nelY , nelZ ,48,44,16,10);
solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, 293.15*ones(size(solidThermoObject.meshObject.nodes,1),1)];

solidThermoObject.elementDisplacementType = 'displacement';
solidThermoObject.materialObject.name = 'Hooke3D';  %   Hooke3D % HeatEquation % Ibrahimbegovic
solidThermoObject.materialObject.rho = 8.92e+3 ;    % density
solidThermoObject.materialObject.E = E;
solidThermoObject.materialObject.mu = E/(2*(1+ny));
solidThermoObject.materialObject.lambda = E*ny/((1+ny)*(1-2*ny));
solidThermoObject.materialObject.kappa = 400;   % thermal conductivity
solidThermoObject.materialObject.alpha = alpha*eye(3,3);    % coupling parameter
solidThermoObject.materialObject.thetaR = [293.15 293.15 293.15 293.15 293.15 293.15 293.15 293.15];         % reference Temperature
solidThermoObject.materialObject.k0 = 385;      % heat capacity
solidThermoObject.dimension = 3;
solidThermoObject.shapeFunctionObject.order = 1;
solidThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
solidThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidThermoObject.numericalTangentObject.showDifferences = false;

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidThermoObject.meshObject.nodes(:,1)==0);
dirichletObject1.nodalDof = 1:3;
dirichletObject1.masterObject = solidThermoObject;
dirichletObject.timeFunction = str2func('@(t) 0');

% mechanical load
neumannObject = neumannClass(dofObject);
neumannObject.masterObject = solidThermoObject;
neumannObject.loadVector = [0;0;0];
neumannObject.loadGeometry = 'area';
neumannObject.timeFunction = str2func('@(t) t*(t<(3600*20))');
neumannObject.meshObject.edof = edofNeumann;

% thermal load
neumannObject2 = neumannClass(dofObject);
neumannObject2.masterObject = solidThermoObject;
neumannObject2.loadGeometry = 'area';
neumannObject2.loadPhysics = 'thermal';
neumannObject2.loadVector = 1000;
neumannObject2.timeFunction = str2func('@(t) (t<(3600*10 - 3600*30/200))' );
neumannObject2.meshObject.edof = edofNeumann;

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
