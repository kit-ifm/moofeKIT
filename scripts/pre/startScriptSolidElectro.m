%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solidElectro';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 30;
setupObject.totalTime = 3;
setupObject.plotObject.flag = false;
setupObject.plotObject.makeMovie = true;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidElectroObject = solidElectroClass(dofObject);
% solidElectroObject = solidClass(dofObject);
[solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof] = meshGeneratorCube(4,4,4,4,4,4,1,false);
solidElectroObject.meshObject.nodes = solidElectroObject.meshObject.nodes + 2;
solidElectroObject.meshObject.nodes = [solidElectroObject.meshObject.nodes, zeros(size(solidElectroObject.meshObject.nodes,1),1)];
% solidElectroObject.elementDisplacementType = 'displacementSC';
% solidElectroObject.materialObject.name = 'NeoHooke';
solidElectroObject.elementDisplacementType = 'mixedD_SC';
solidElectroObject.materialObject.name = 'MooneyRivlin';
solidElectroObject.materialObject.rho = 0.1;
% Material Parameters
solidElectroObject.materialObject.a = 1e3;
solidElectroObject.materialObject.b = 1e3;
solidElectroObject.materialObject.c = 1e2;
solidElectroObject.materialObject.d = 2*(solidElectroObject.materialObject.a + 2*solidElectroObject.materialObject.b);
solidElectroObject.materialObject.epsilon0 = 0;
solidElectroObject.materialObject.e0 = 8.854*10^(-12);
solidElectroObject.materialObject.e1 = 4;
solidElectroObject.materialObject.e2 = 0;
solidElectroObject.materialObject.ee = 0;
solidElectroObject.materialObject.rhoSource = 0;
solidElectroObject.materialObject.timeFunctionRhoSource = @(t) 0;
% 
solidElectroObject.dimension = 3;
solidElectroObject.shapeFunctionObject.order = 1;
solidElectroObject.shapeFunctionObject.numberOfGausspoints = 8;
solidElectroObject.mixedFEObject.shapeFunctionObject.order = solidElectroObject.shapeFunctionObject.order - 1;
solidElectroObject.mixedFEObject.condensation = false;

% dirichletObject1 = dirichletClass(dofObject);
% dirichletObject1.nodeList = find(solidElectroObject.nodes(:,3) == 0);
% dirichletObject1.nodalDof = 1:3;
% dirichletObject1.masterObject = solidElectroObject;
% dirichletObject1.timeFunction = str2func('@(t,XYZ) 0');

% dirichletElectro1 = dirichletClass(dofObject);
% dirichletElectro1.masterObject = solidElectroObject;
% dirichletElectro1.nodeList = find(solidElectroObject.nodes(:,3)==0);
% dirichletElectro1.nodalDof = 4;
% dirichletElectro1.timeFunction = str2func('@(t,addfield1) 0');

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidElectroObject.meshObject.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = solidElectroObject;
boundary1.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

% dirichletElectro3 = dirichletClass(dofObject);
% dirichletElectro3.masterObject = solidElectroObject;
% dirichletElectro3.nodeList = find(solidElectroObject.nodes(:,1));
% dirichletElectro3.nodalDof = 4;
% dirichletElectro3.timeFunction = str2func('@(t,phi) 0');

dirichletElectro3 = dirichletClass(dofObject);
dirichletElectro3.masterObject = solidElectroObject;
% dirichletElectro3.nodeList = find((solidElectroObject.nodes(:,1)==4)&(solidElectroObject.nodes(:,3)>=2));
dirichletElectro3.nodeList = find((solidElectroObject.meshObject.nodes(:,1)==4));
dirichletElectro3.nodalDof = 4;
dirichletElectro3.timeFunction = str2func('@(t,phi) t*25*(t<1) + 25*(t>=1)');

dirichletElectro4 = dirichletClass(dofObject);
dirichletElectro4.masterObject = solidElectroObject;
% dirichletElectro4.nodeList = find((solidElectroObject.nodes(:,1)==0)&(solidElectroObject.nodes(:,3)>=2));
dirichletElectro4.nodeList = find((solidElectroObject.meshObject.nodes(:,1)==0));
dirichletElectro4.nodalDof = 4;
dirichletElectro4.timeFunction = str2func('@(t,phi) -t*25*(t<1) - 25*(t>=1)');

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
omegaEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'omegaEnergy');
ePotMechanical = getEnergy(dofObject.postDataObject,dofObject,setupObject,'ePotMechanical');
ePotElectrical = getEnergy(dofObject.postDataObject,dofObject,setupObject,'ePotElectrical');
ePotD0timesGRADphi = getEnergy(dofObject.postDataObject,dofObject,setupObject,'ePotD0timesGRADphi');
figure; 
plot(timeVector,kineticEnergy + ePotMechanical + ePotElectrical + ePotD0timesGRADphi);
figure;
plot(solidElectroObject)
% filename = 'execute';
% VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))

