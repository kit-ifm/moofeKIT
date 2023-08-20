% MICROFLUIDICPUMPINGDEVICE Script for preprocessing a dynamic coupled 
% electro-thermo-mechanical simulation.
% 
% FORMULATION
% Mixed finite element formulation for nonlinear electro-thermo-mechanical 
% processes is pursued. This is accomplished under the assumption of a 
% hyperelastic, isotropic material and based on a polyconvex inspired 
% internal energy function. The discretization in time is pursued with an 
% energy and momentum conserving scheme. 
% 
% 
% REFERENCE
% paper in preparation
% https://doi.org/10.1016/j.cma.2021.114298
% https://doi.org/10.1007/BF00913408
% 
% SEE ALSO
% LShape, 
% cooksMembrane
% 
% CREATOR(S) 
% Marlon Franke
 
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = mfilename;
setupObject.saveObject.saveData = true;
setupObject.totalTimeSteps = 1000;
setupObject.saveObject.saveForTimeSteps = [1 2:2:setupObject.totalTimeSteps];
% setupObject.saveObject.saveForTimeSteps = 1:setupObject.totalTimeSteps;
setupObject.totalTime = 100;
setupObject.newton.tolerance = 1e-6;
setupObject.newton.maximumSteps = 10;
setupObject.newton.enforceIteration = true;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'stress'; %'temp'; 'phi'; 'D1';
setupObject.plotObject.view = [-0.4,90];
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

%% dofObject required object for dof and object handling
dofObject = dofClass;   

%% mesh
if 0 % H1H0
%     abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH1.inp');
    abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH1Coarse2.inp');
    shapeFunctionOrder = 1;
%     mixedFEShapeFunctionOrder = 0;
    mixedFEShapeFunctionOrder = 1;
    numberOfGausspoints = 8;
else % 
%     abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH2Serendipity.inp');
%     abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH2SerendipityMedium.inp');
    abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH2SerendipityCoarse.inp');
%     abaqusMeshData = abaqusInputFileConverter('microFluidicPumpingDeviceH2SerendipityCoarse2.inp');
    shapeFunctionOrder = 2;
    mixedFEShapeFunctionOrder = 1;
    numberOfGausspoints = 27;
end

solidElectroThermoObject = solidElectroThermoClass(dofObject);
solidElectroThermoObject.meshObject.nodes = abaqusMeshData.qR;
solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes*1e-3;
nnodes=size(solidElectroThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293*ones(nnodes,1);
initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
solidElectroThermoObject.meshObject.edof = abaqusMeshData.edof;
solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';
solidElectroThermoObject.elementDisplacementType = 'mixedSC';
solidElectroThermoObject.materialObject.rho = 1000*1e-9;                                          % mass-density
solidElectroThermoObject.materialObject.a = 25000*1e-3;                                         % a
solidElectroThermoObject.materialObject.b = 50000*1e-3;                                         % b
solidElectroThermoObject.materialObject.c1 = 500000*1e-3;                                         % c
solidElectroThermoObject.materialObject.c2 = 5209*1e-3;                                          % e
solidElectroThermoObject.materialObject.d1 = 2*(solidElectroThermoObject.materialObject.a + 2*solidElectroThermoObject.materialObject.b);   % d
solidElectroThermoObject.materialObject.d2 = 0;
solidElectroThermoObject.materialObject.e0 = 8.8541*10^(-12)*1e-9;                             	% epsilon_0
solidElectroThermoObject.materialObject.e1 = 4;                                             % epsilon_r
solidElectroThermoObject.materialObject.e2 = 0;
solidElectroThermoObject.materialObject.ee = 0;
solidElectroThermoObject.materialObject.kappa = 1500*1e-3;                                        % heat capacity
solidElectroThermoObject.materialObject.strongcomp = false;
solidElectroThermoObject.materialObject.beta = 2.233*10^(-4);                             % coupling parameter
solidElectroThermoObject.materialObject.thetaR = 293.15;                                  % Reference Temperature
solidElectroThermoObject.materialObject.k0 = 0.23*1e3;                                          % thermal conductivity
solidElectroThermoObject.materialObject.wk = 0;
solidElectroThermoObject.materialObject.rhoSource = 0;
solidElectroThermoObject.materialObject.RSource = 0;
solidElectroThermoObject.materialObject.timeFunctionRhoSource = @(t) 0;
solidElectroThermoObject.dimension = 3;
solidElectroThermoObject.shapeFunctionObject.order = shapeFunctionOrder;
solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints;
solidElectroThermoObject.mixedFEObject.condensation = true;
solidElectroThermoObject.mixedFEObject.typeShapeFunction = 'sameOrder';
solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = mixedFEShapeFunctionOrder;
solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidElectroThermoObject.numericalTangentObject.showDifferences = false;

% electrical dirichlet boundary condition, DirichletElectricalLower
dirichletObject1a = dirichletClass(dofObject);
dirichletObject1a.nodeList = abaqusMeshData.subsets(1).edof;
dirichletObject1a.nodalDof = 4;
dirichletObject1a.masterObject = solidElectroThermoObject;
phi0 = 7.5e8;
dirichletObject1a.timeFunction = @(t,Z) microFluidicPumpingDeviceMyFunction(t);
% dirichletObject1a.timeFunction = @(t,Z) (phi0*sin(pi/2*t/bcTimeEnd(2))).*(t >= bcTimeEnd(1)).*(t <= bcTimeEnd(2)) + ...
%                                         phi0*(t > bcTimeEnd(2)).*(t <= bcTimeEnd(3)) + ...
%                                        (phi0*cos(pi/2*(t-bcTimeEnd(3))/(bcTimeEnd(4)-bcTimeEnd(3)))).*(t > bcTimeEnd(3)).*(t <= bcTimeEnd(4)) + ...
%                                        (phi0*sin(pi/2*(t-bcTimeEnd(9))/(bcTimeEnd(10)-bcTimeEnd(9)))).*(t > bcTimeEnd(9)).*(t <= bcTimeEnd(10)) + ... 
%                                        phi0*(t > bcTimeEnd(10)).*(t <= bcTimeEnd(11));
% nodes3 = nodes(dirichletObject1a.nodeList,1:3); plot3(nodes3(:,1),nodes3(:,2),nodes3(:,3),'*')

% electrical dirichlet boundary condition, DirichletElectricalLower
dirichletObject1b = dirichletClass(dofObject);
dirichletObject1b.nodeList = abaqusMeshData.subsets(3).edof;
dirichletObject1b.nodalDof = 4;
dirichletObject1b.masterObject = solidElectroThermoObject;
dirichletObject1b.timeFunction = @(t,Z) -microFluidicPumpingDeviceMyFunction(t+16);
% dirichletObject1b.timeFunction = @(t,Z) (0*sin(pi/2*t/bcTimeEnd(2))).*(t >= 0).*(t <= bcTimeEnd(2)) + 0*(t > bcTimeEnd(2)).*(t <= bcTimeEnd(3)) + ...
%                                        (0*cos(pi/2*(t-bcTimeEnd(3))/(bcTimeEnd(4)-bcTimeEnd(3)))).*(t > bcTimeEnd(3)).*(t <= bcTimeEnd(4)) + 0*(t > bcTimeEnd(4)).*(t <= bcTimeEnd(5)) + ...
%                                        (-phi0*sin(pi/2*(t-bcTimeEnd(5))/(bcTimeEnd(6)-bcTimeEnd(5)))).*(t > bcTimeEnd(5)).*(t <= bcTimeEnd(6)) + -phi0*(t > bcTimeEnd(6)).*(t <= bcTimeEnd(7)) + ...
%                                        (-phi0*cos(pi/2*(t-bcTimeEnd(7))/(bcTimeEnd(8)-bcTimeEnd(7)))).*(t > bcTimeEnd(7)).*(t <= bcTimeEnd(8)) + 0*(t > bcTimeEnd(8)).*(t <= bcTimeEnd(9)) + ...
%                                        (0*sin(pi/2*(t-bcTimeEnd(9))/(bcTimeEnd(10)-bcTimeEnd(9)))).*(t > bcTimeEnd(9)).*(t <= bcTimeEnd(10)) + 0*(t > bcTimeEnd(10)).*(t <= bcTimeEnd(11));
% nodes3 = nodes(dirichletObject1b.nodeList,1:3); plot3(nodes3(:,1),nodes3(:,2),nodes3(:,3),'*')

% electrical dirichlet boundary condition, DirichletElectricalMid
dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = abaqusMeshData.subsets(2).edof;
dirichletObject2.nodalDof = 4;
dirichletObject2.masterObject = solidElectroThermoObject;
% nodes3 = nodes(dirichletObject2.nodeList,1:3); plot3(nodes3(:,1),nodes3(:,2),nodes3(:,3),'*')

% mechanical dirichlet boundary condition, DirichletFixed
dirichletObject4 = dirichletClass(dofObject);
dirichletObject4.nodeList = abaqusMeshData.subsets(4).edof;
dirichletObject4.nodalDof = 1:3;
dirichletObject4.masterObject = solidElectroThermoObject;
% nodes4 = nodes(dirichletObject4.nodeList,1:3); plot3(nodes4(:,1),nodes4(:,2),nodes4(:,3),'*')

% mechanical dirichlet boundary condition, DirichletSymmetryX
dirichletObject5 = dirichletClass(dofObject);
dirichletObject5.nodeList = abaqusMeshData.subsets(5).edof;
dirichletObject5.nodalDof = 1;
dirichletObject5.masterObject = solidElectroThermoObject;
% nodes5 = nodes(dirichletObject5.nodeList,1:3); plot3(nodes5(:,1),nodes5(:,2),nodes5(:,3),'*')

% mechanical dirichlet boundary condition, DirichletSymmetryY
dirichletObject6 = dirichletClass(dofObject);
dirichletObject6.nodeList = abaqusMeshData.subsets(6).edof;
dirichletObject6.nodalDof = 2;
dirichletObject6.masterObject = solidElectroThermoObject;
% nodes6 = nodes(dirichletObject6.nodeList,1:3); plot3(nodes6(:,1),nodes6(:,2),nodes6(:,3),'*')

% mechanical dirichlet boundary condition, DirichletSymmetryZ
dirichletObject7 = dirichletClass(dofObject);
dirichletObject7.nodeList = abaqusMeshData.subsets(7).edof;
dirichletObject7.nodalDof = 3;
dirichletObject7.masterObject = solidElectroThermoObject;
% nodes7 = nodes(dirichletObject7.nodeList,1:3); plot3(nodes7(:,1),nodes7(:,2),nodes7(:,3),'*')

% thermal neumann boundary condition, NeumannThermalGeometry
neumannObject1 = neumannClass(dofObject);
neumannObject1.field = 'thermal';
neumannObject1.masterObject = solidElectroThermoObject;
neumannObject1.energy = 0.8;
neumannObject1.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);
bcTimeEnd = 0:4:100;
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd(2)/2)+(bcTimeEnd(2)-t).*(t > bcTimeEnd(2)/2 & t <= bcTimeEnd(2));
neumannObject1.meshObject.edof = abaqusMeshData.subsets(23).edof;
% nodes8 = nodes(unique(neumannObject1.meshObject.edof(:)),1:3); plot3(nodes8(:,1),nodes8(:,2),nodes8(:,3),'*')


%% Test thermal dirichlet boundary condition on all nodes -> does not work properly!
% dirichletObject9 = dirichletClass(dofObject);
% nodeList2 = 1:size(solidElectroThermoObject.meshObject.nodes,1);
% dirichletObject9.nodeList = nodeList2;
% dirichletObject9.nodalDof = 5;
% dirichletObject9.masterObject = solidElectroThermoObject;
% dirichletObject9.timeFunction = str2func('@(t,theta) 293.15');  

% % electrical neumann boundary condition, DirichletElectricalLower
% neumannObject0 = neumannClass(dofObject);
% neumannObject0.field = 'electrical';
% neumannObject0.omega_0 = 0.01;
% neumannObject0.meshObject.edof = abaqusMeshData.subsets(20).edof;
% neumannObject0.masterObject = solidElectroThermoObject;
% neumannObject0.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
% neumannObject0.shapeFunctionObject.order = solidElectroThermoObject.shapeFunctionObject.order;
% neumannObject0.shapeFunctionObject.numberOfGausspoints = 2^(solidElectroThermoObject.dimension-1);

% % zero dirichlet boundary condition of all non-eap-nodes
% dirichletObject3 = dirichletClass(dofObject);
% nodeList = 1:size(solidElectroThermoObject.meshObject.nodes,1);
% nodes = solidElectroThermoObject.meshObject.nodes;
% % nodeListZ = (nodes(:,3) <= 99*1e-6).*(nodes(:,3) > 81*1e-6);
% % nodeListXY = sqrt(nodes(:,1).^2+nodes(:,2).^2) < 460.5*1e-6;
% nodeListZ = (nodes(:,3) <= 99).*(nodes(:,3) > 81);
% nodeListXY = sqrt(nodes(:,1).^2+nodes(:,2).^2) < 460.5;
% nodeListXYZ = find(nodeListZ.*nodeListXY);
% nonNodeListTemp = unique([abaqusMeshData.subsets(1).edof; nodeListXYZ]);
% nodeList(nonNodeListTemp) = [];
% % assert(numel(nodeListXYZ)==numel(unique(abaqusMeshData.subsets(1).edof(:))),'number of nodes are not consistent')
% dirichletObject3.nodeList = nodeList;
% dirichletObject3.nodalDof = 4;
% dirichletObject3.masterObject = solidElectroThermoObject;
% % nodes3 = nodes(nodeList,1:3); plot3(nodes3(:,1),nodes3(:,2),nodes3(:,3),'*')
% % nodes3b = nodes(nonNodeListTemp,1:3); plot3(nodes3b(:,1),nodes3b(:,2),nodes3(:,3),'*')

% % thermal dirichlet boundary condition, DirichletThermal
% dirichletObject8 = dirichletClass(dofObject);
% dirichletObject8.nodeList = abaqusMeshData.subsets(7).edof;
% dirichletObject8.nodalDof = 5;
% dirichletObject8.masterObject = solidElectroThermoObject;
% dirichletObject8.timeFunction = str2func('@(t,theta) 293.15+50');  