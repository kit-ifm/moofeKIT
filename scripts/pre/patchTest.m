%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solidElectroThermoPatchTestDistorted';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.newton.tolerance = 1e-6;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress.component = -1;
setupObject.plotObject.makeMovie = false;
% setupObject.integrator = 'DiscreteGradient';
setupObject.integrator = 'Endpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% solidElectroThermoObject = solidElectroThermoClass(dofObject);
solidElectroThermoObject = solidClass(dofObject);
[solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1,1,1,2,2,2,2,true);
solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
% solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
% numberOfElements = 2;
% [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof,edofneumannTet] = triquadraticTetrahedralBrick(numberOfElements, numberOfElements, numberOfElements, 1, 1, 1, 1);
% solidElectroThermoObject.meshObject.nodes(:,3) = solidElectroThermoObject.meshObject.nodes(:,3) + 1;
% solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 11;

% solidElectroThermoObject.meshObject.nodes = [0 0 0;                  %1
%                             1 0 0;                  %2
%                             1 1 0;                  %3
%                             0 1 0;                  %4
%                             0 0 1;                  %5
%                             1 0 1;                  %6
%                             1 1 1;                  %7
%                             0 1 1;                  %8
%                             0.249 0.342 0.192;      %9
%                             0.826 0.288 0.288;      %10
%                             0.850 0.649 0.263;      %11
%                             0.273 0.750 0.230;      %12
%                             0.320 0.186 0.643;      %13
%                             0.677 0.305 0.683;      %14
%                             0.788 0.693 0.644;      %15
%                             0.165 0.745 0.702;      %16
%                             0.5 0 0;                %17
%                             1 0.5 0;                %18
%                             0.5 1 0;                %19
%                             0 0.5 0;                %20
%                             0.5 0 1;                %21
%                             1 0.5 1;                %22
%                             0.5 1 1;                %23
%                             0 0.5 1;                %24
%                             0 0 0.5;                %25
%                             1 0 0.5;                %26
%                             1 1 0.5;                %27
%                             0 1 0.5;                %28
%                             0.16 0.093 0.8215;       %29
%                             0.8385    0.1525    0.8415;
%                             0.9130    0.1440    0.1440;
%                             0.1245    0.1710    0.0960;
%                             0.7515    0.2965    0.4855;
%                             0.5375    0.3150    0.2400;
%                             0.2845    0.2640    0.4175;
%                             0.4985    0.2455    0.6630;
%                             0.8380    0.4685    0.2755;
%                             0.2610    0.5460    0.2110;
%                             0.2425    0.4655    0.6725;
%                             0.7325    0.4990    0.6635;      %40
%                             0.8190    0.6710    0.4535;
%                             0.5615    0.6995    0.2465;
%                             0.2190    0.7475    0.4660;
%                             0.4765    0.7190    0.6730;
%                             0.8940    0.8465    0.8220;
%                             0.9250    0.8245    0.1315;
%                             0.1365    0.8750    0.1150;
%                             0.0825    0.8725    0.8510];      %48
% solidElectroThermoObject.meshObject.edof =[
%                             1 2 3 4 9 10 11 12 17 18 19 20 34 37 42 38 32 31 46 47;
%                             16 15 11 12 8 7 3 4 44 41 42 43 23 27 19 28 48 45 46 47;
%                             16 8 5 13 15 7 6 14 48 24 29 39 45 22 30 40 44 23 21 36;
%                             9 13 5 1 10 14 6 2 35 29 25 32 33 30 26 31 34 36 21 17;
%                             11 15 14 10 3 7 6 2 41 40 33 37 27 22 26 18 46 45 30 31;
%                             4 8 5 1 12 16 13 9 28 24 25 20 43 39 35 38 47 48 29 32;
%                             12 16 13 9 11 15 14 10 43 39 35 38 41 40 33 37 42 44 36 34];
nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);
initialElectricalPotentialField = zeros(nnodes,1);
% solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.name = 'MooneyRivlin';
solidElectroThermoObject.materialObject.rho = 0;                                            % mass-density
solidElectroThermoObject.materialObject.a = 831.25;                                         % a
solidElectroThermoObject.materialObject.b = 166.25;                                         % b
solidElectroThermoObject.materialObject.c = 10000;                                         % c
solidElectroThermoObject.materialObject.c2 = 5209;                                          % e
solidElectroThermoObject.materialObject.d = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12);                                % epsilon_0
solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 100;                                        % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                               % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 293.15;                                    % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 10;                                            % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = 2;
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
% strainEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
% externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
% %figure;
% %plot(timeVector,kineticEnergy + strainEnergy);
% %plot(solidThermoObject)
% % filename = 'execute';
% % VTKPlot(filename,'unstructured_grid',part.qN1(:,1:3),part.edof,'scalars','temperature',part.qN1(:,4))
