%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solidElectroThermoPatchTestDistorted';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'phi';
setupObject.plotObject.stress.component = -1;
setupObject.plotObject.makeMovie = false;
setupObject.integrator = 'DiscreteGradient';
% setupObject.integrator = 'Endpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidElectroThermoObject = solidElectroThermoClass(dofObject);
% [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1,1,1,2,2,2,2,true);
% solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
[solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshPatchTestDistorted(1,1,1,2,true);

nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);
initialElectricalPotentialField = zeros(nnodes,1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
solidElectroThermoObject.materialObject.rho = 0;                                            % mass-density
solidElectroThermoObject.materialObject.a = 25000;                                         % a
solidElectroThermoObject.materialObject.b = 50000;                                         % b
solidElectroThermoObject.materialObject.c1 = 500000;                                         % c
solidElectroThermoObject.materialObject.c2 = 5209;                                          % e
solidElectroThermoObject.materialObject.d1 = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12);                                % epsilon_0
solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 1500;                                        % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                               % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 293.15;                                    % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 0.23;                                            % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
solidElectroThermoObject.mixedFEObject.condensation = false;
solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = 1;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = false;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = solidElectroThermoObject;

boundary1(2) = dirichletClass(dofObject);
boundary1(2).nodeList = find(solidElectroThermoObject.meshObject.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = solidElectroThermoObject;

boundary1(3) = dirichletClass(dofObject);
boundary1(3).nodeList = find(solidElectroThermoObject.meshObject.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = solidElectroThermoObject;

boundaryK2 = dirichletClass(dofObject);
boundaryK2.nodeList = find(solidElectroThermoObject.meshObject.nodes(:,3) == 1);
boundaryK2.nodalDof = 3;
boundaryK2.masterObject = solidElectroThermoObject;
boundaryK2.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');  

%% solver
tic
dofObject = runNewton(setupObject,dofObject);
toc

%% postprocessing - energy
% f = figure;
% plot(solidElectroThermoObject, setupObject);
% xlim([0 inf]);
% ylim([0 inf]);
% zlim([0 inf]);
% xticks([]);
% yticks([]);
% zticks([]);
% f.Position = [100 100 320 210];
% % zlim([0 1])
% % xlabel('x');
% % ylabel('y');
% % zlabel('z');
% caxis([-1 1]);
% campos([-5,-6,5]);
% %set(gcf, 'Renderer', 'painter');
% % export_fig('test','-eps');
% 
% timeVector = getTime(dofObject.postDataObject,setupObject);
% kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
% strainEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
% externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
% %figure; 
% %plot(timeVector,kineticEnergy + strainEnergy);
% %plot(solidThermoObject)
% % filename = 'execute';
% % VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))
