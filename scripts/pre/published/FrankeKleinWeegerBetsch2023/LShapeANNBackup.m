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

run('../../startUpMoofeKIT.m')

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShape';
setupObject.saveObject.saveData = true;
setupObject.totalTimeSteps = 160;
% setupObject.totalTimeSteps = 200;
setupObject.totalTime = 100;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-6;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
mesh = 'H1H0'; 
% mesh = 'H2H1';
solidObject = solidClass(dofObject);
if strcmpi(mesh,'H1H0')
%     abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
%     abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
elseif strcmpi(mesh,'H2H1')
    abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
end
solidObject.meshObject.nodes = abaqusMeshData.qR;
solidObject.meshObject.edof = abaqusMeshData.edof;
% 
solidObject.materialObject.name = 'ANNMooneyRivlin';    % IC, IIC, c
% solidObject.materialObject.name = 'ANNMooneyRivlin2';    % C, G, c %TODO: implement analytical tangent
% solidObject.materialObject.name = 'ANNMooneyRivlin3';    % IC, IIC, J, -J
% solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
solidObject.elementDisplacementType = 'displacementSC';
%     solidObject.elementDisplacementType = 'mixedSC';
% 
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = false;
% 
bulk    = 5209;                                                 % K
shear   = 997.5;                                                % G
mu      = shear;                                                % second lame parameter
c1      = 5/6*mu;                                               % a
c2      = 1/6*mu;                                               % b
c       = 0; 
solidObject.materialObject.a = c1/2;
solidObject.materialObject.b = c2/2;
solidObject.materialObject.c = c;
solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
solidObject.mixedFEObject.condensation = true;
solidObject.materialObject.rho = 100;
solidObject.dimension = 3;
% solidObject.shapeFunctionObject.numberOfGausspoints = 27;

%% ANN
% already trained ANN model
solidObject.artificialNeuralNetworkObject.computeActivationFunction('softPlus');
% solidObject.artificialNeuralNetworkObject.alpha = 5e2;
neuronName = '[8]/'; % 8, 16, 32, 64, 128
% ANNlocationFolder = '~/Projekte/moofeKIT/continuumClasses/@solidClass/ANNData/normalized/';
ANNlocationFolder = '~/Projekte/moofeKIT/continuumClasses/@solidClass/ANNData/normalized2/';
solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = strcat(ANNlocationFolder,neuronName);
% ANNlocationFolder = '~/Projekte/moofeKIT/continuumClasses/@solidClass/ANNData/unnormalized/';
% solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = ANNlocationFolder;
solidObject.artificialNeuralNetworkObject.readData.biasNameCell = {'b1','b2'};
solidObject.artificialNeuralNetworkObject.readData.weightNameCell = {'w1','w2'};
solidObject.artificialNeuralNetworkObject.read('bias');
solidObject.artificialNeuralNetworkObject.read('weight');
solidObject.artificialNeuralNetworkObject.normalizationFactor = readmatrix(strcat(ANNlocationFolder,'normalizationFactor.txt'));

% bcTimeEnd = 1;
FA = [256;512;768]*4;
bcTimeEnd = 5;
% FA = [256;512;768];
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.typeOfLoad = 'deadLoad';
neumannObject1.masterObject = solidObject;
neumannObject1.forceVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
% TODO scale by bcTimeEnd/2
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.typeOfLoad = 'deadLoad';
neumannObject2.masterObject = solidObject;
neumannObject2.forceVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction =  @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
% TODO scale by bcTimeEnd/2
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

%% solver

%try
    dofObject = runNewton(setupObject,dofObject);
%end

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
try
    [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
    [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
end
strainEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = strainEnergy(1:length(kineticEnergy)) + kineticEnergy;
if strcmpi(setupObject.integrator,'DiscreteGradient')
    tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
    totalEnergyDiff = totalEnergy(tStartDiff+1:end) - totalEnergy(tStartDiff:end-1);
    figure;
    plot(timeVector(tStartDiff+1:end), totalEnergyDiff);
    matlab2tikz(['diffEnergy',solidObject.materialObject.name,setupObject.integrator,'.tikz'],'width','\figW','height','\figH')
end
try
%     figure;
%     plot(timeVector,linearMomentum);
%     figure;
%     plot(timeVector,totalLinearMomentum);
    figure;
    plot(timeVector,angularMomentum);
    xlabel t; ylabel J; legend('J_x','J_y','J_z'); set(gca,'fontsize',20)
    matlab2tikz(['angularMomentum',solidObject.materialObject.name,setupObject.integrator,'.tikz'],'width','\figW','height','\figH')
    figure;
    plot(timeVector,totalAngularMomentum);
    xlabel t; ylabel J; legend('EM'); set(gca,'fontsize',20)
    matlab2tikz(['totalAngularMomentum',solidObject.materialObject.name,setupObject.integrator,'.tikz'],'width','\figW','height','\figH')
end
figure;
plot(timeVector, totalEnergy);
xlabel t; ylabel E; set(gca,'fontsize',20)
matlab2tikz(['energy',solidObject.materialObject.name,setupObject.integrator,'.tikz'],'width','\figW','height','\figH')
% 
figure;
plot(timeVector, setupObject.newton.step','x');
xlabel t; ylabel Iterations; set(gca,'fontsize',20)
matlab2tikz(['iterations',solidObject.materialObject.name,setupObject.integrator,'.tikz'],'width','\figW','height','\figH')
