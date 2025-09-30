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

clear all;

% run('../../startUpMoofeKIT.m')

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShape';
setupObject.saveObject.saveData = false;
% setupObject.totalTimeSteps = 125;
% setupObject.totalTimeSteps = 320;% setupObject.totalTime = 200;
% setupObject.totalTimeSteps = 250;
setupObject.totalTimeSteps = 200;
% setupObject.totalTime = 1;
setupObject.totalTime = 200;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.view = [-0.4,90];
setupObject.plotObject.view = [0,90];
setupObject.plotObject.colorBarLimits = [0 600];
setupObject.plotObject.border = 0.4;
setupObject.plotObject.lineWidth = 2;
setupObject.plotObject.savePlot.flag = true;
setupObject.plotObject.savePlot.name = 'LShapeGT';
setupObject.plotObject.savePlot.type = '-dpng';
setupObject.newton.tolerance = 5e-4;
% setupObject.newton.tolerance = 1e-5; % not good for ANN model
setupObject.newton.tolerance = 1e-5;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';
%
neuronNameVector = [4,8,16,32,64,128];
neuronNumberCounter = 2;
dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
% mesh = 'H1H0';
mesh = 'H1';
% mesh = 'H2H1';
solidObject = solidClass(dofObject);
if strcmpi(mesh,'H1')
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.elementDisplacementType = 'displacementSC';
    %     solidObject.elementNameAdditionalSpecification = 'IIIJ';
    % solidObject.numericalTangentObject.computeNumericalTangent = true;
    % solidObject.numericalTangentObject.type = 'complex';
elseif strcmpi(mesh,'H1H0')
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.elementDisplacementType = 'mixedSC';
elseif strcmpi(mesh,'H2H1')
    abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.mixedFEObject.condensation = true;
    solidObject.elementDisplacementType = 'mixedSC';
end
solidObject.meshObject.nodes = abaqusMeshData.qR;
solidObject.meshObject.edof = abaqusMeshData.edof;

solidObject.materialObject.name = 'ANN';
% solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
%
bulk    = 5209;                                                 % K
shear   = 997.5;                                                % G
mu      = shear;                                                % second lame parameter
a      = 5/6*mu;                                               % a
b      = 1/6*mu;                                               % b
c       = 10000;
solidObject.materialObject.a = a;
solidObject.materialObject.b = b;
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
neuronName = strcat('[',num2str(neuronNameVector(neuronNumberCounter)),']/');
% ANNlocationFolder = 'externalData/ANN/solid/normalized/';
% ANNlocationFolder = 'externalData/ANN/solid/normalized2/';
ANNlocationFolder = 'externalData/ANN/solid/normalized3/';
solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = strcat(ANNlocationFolder,neuronName);
% ANNlocationFolder = 'externalData/ANN/solid/unnormalized/';
% solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = ANNlocationFolder;
solidObject.artificialNeuralNetworkObject.readData.biasNameCell = {'b1','b2'};
solidObject.artificialNeuralNetworkObject.readData.weightNameCell = {'w1','w2'};
solidObject.artificialNeuralNetworkObject.read('bias');
solidObject.artificialNeuralNetworkObject.read('weight');
solidObject.artificialNeuralNetworkObject.normalizationFactor = readmatrix(strcat(ANNlocationFolder,'normalizationFactor.txt'));
solidObject.artificialNeuralNetworkObject.A = readmatrix(strcat(ANNlocationFolder,neuronName,'A.txt'));

FA = [256;512;768];
bcTimeEnd = 5;
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.loadType = 'deadLoad';
neumannObject1.loadGeometry = 'area';
neumannObject1.masterObject = solidObject;
neumannObject1.loadVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
% TODO scale by bcTimeEnd/2
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.loadType = 'deadLoad';
neumannObject2.loadGeometry = 'area';
neumannObject2.masterObject = solidObject;
neumannObject2.loadVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction =  @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
% TODO scale by bcTimeEnd/2
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

%% solver
% try
dofObject = runNewton(setupObject,dofObject);
% end

%% postprocessing - energy
nameFile = [solidObject.materialObject.name,setupObject.integrator,'DT',num2str(setupObject.timeStepSize)];
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
try
    [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
    [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
end
strainEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = strainEnergy(1:length(kineticEnergy)) + kineticEnergy + externalEnergy(1:length(kineticEnergy));
if strcmpi(setupObject.integrator,'DiscreteGradient')
    tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
    totalEnergyDiff = totalEnergy(tStartDiff+2:end) - totalEnergy(tStartDiff+1:end-1);
    figure;
    plot(timeVector(tStartDiff+1:end-1), totalEnergyDiff);
    matlab2tikz(['diffEnergy',nameFile,'.tikz'],'width','\figW','height','\figH')
end
try
    %     figure;
    %     plot(timeVector,linearMomentum);
    %     figure;
    %     plot(timeVector,totalLinearMomentum);
    figure;
    plot(timeVector,angularMomentum);
    xlabel t; ylabel J; legend('J_x','J_y','J_z'); set(gca,'fontsize',20)
    matlab2tikz(['angularMomentum',nameFile,'.tikz'],'width','\figW','height','\figH')
    figure;
    plot(timeVector,totalAngularMomentum);
    xlabel t; ylabel J; legend('EM'); set(gca,'fontsize',20)
    matlab2tikz(['totalAngularMomentum',nameFile,'.tikz'],'width','\figW','height','\figH')
end
figure;
plot(timeVector, totalEnergy);
xlabel t; ylabel E; set(gca,'fontsize',20)
matlab2tikz(['energy',nameFile,'.tikz'],'width','\figW','height','\figH')
%
figure;
plot(timeVector(2:end), setupObject.newton.step','ro');
xlabel t; ylabel Iterations; set(gca,'fontsize',20)
matlab2tikz(['iterations',nameFile,'.tikz'],'width','\figW','height','\figH')
