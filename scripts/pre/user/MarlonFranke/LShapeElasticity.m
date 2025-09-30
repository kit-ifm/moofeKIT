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
setupObject.totalTimeSteps = 10000;
% setupObject.totalTimeSteps = 320;% setupObject.totalTime = 200;
% setupObject.totalTimeSteps = 250;
% setupObject.totalTimeSteps = 100;
% setupObject.totalTime = 1;
setupObject.totalTime = 100;
setupObject.plotObject.flag = true;
% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.view = [-0.4,90];
setupObject.plotObject.view = [0,90];
setupObject.plotObject.colorBarLimits = [0 600];
setupObject.plotObject.border = 0.4;
setupObject.plotObject.lineWidth = 2;
setupObject.plotObject.colorscheme = 'turbo';
setupObject.plotObject.steps = 1;
setupObject.plotObject.savePlot.flag = false;
setupObject.plotObject.savePlot.name = 'LShape';
setupObject.plotObject.savePlot.flagPlotColorbar = false;
% setupObject.plotObject.xminxmaxset = [-4 20 -4 14 -4 4];
setupObject.plotObject.savePlot.type = '-dpng';
% setupObject.plotObject.savePlot.type = '-depsc';
setupObject.newton.tolerance = 5e-4;
% setupObject.newton.tolerance = 1e-6; % not good for ANN model<
% setupObject.newton.tolerance = 1e-8;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';
setupObject.integrator = 'LinearImplicit';
%
neuronNameVector = [4,8,16,32,64,128];
neuronNumberCounter = 2;
dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
mesh = 'H1H0';
% mesh = 'H1';
% mesh = 'H2';
% mesh = 'H2H1';
% mesh = 'H2H0';
% mesh = 'H1H1';
% mesh = 'H2H2';
if strcmpi(setupObject.integrator,'LinearImplicit')
    solidObject = solidVelocityClass(dofObject);
else
    solidObject = solidClass(dofObject);
end
if strcmpi(mesh,'H1')
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
    %     abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.elementDisplacementType = 'displacementSC';
    % solidObject.elementNameAdditionalSpecification = 'pH';
    %     solidObject.elementNameAdditionalSpecification = 'IIIJ';
    % solidObject.numericalTangentObject.computeNumericalTangent = true;
    % solidObject.numericalTangentObject.type = 'complex';
elseif strcmpi(mesh,'H2')
    % abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityVeryCoarse.inp');
    % solidObject.elementDisplacementType = 'mixedSC';
    solidObject.elementDisplacementType = 'displacementSC';
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
elseif strcmpi(mesh,'H1H0')
    abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
    % abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    % % % % 
    % abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    % solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.mixedFEObject.typeShapeFunctionData = 0;
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.flagHistoryFields = true;
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    % solidObject.elementNameAdditionalSpecification = 'pHCGc';
    % solidObject.elementNameAdditionalSpecification = 'pHCTilde';
    % solidObject.elementNameAdditionalSpecification = 'pHC';
    solidObject.mixedFEObject.condensation = false;
    % solidObject.mixedFEObject.condensation = true;
elseif strcmpi(mesh,'H2H1')
    % abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
    abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityVeryCoarse.inp');
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.flagHistoryFields = true;
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.condensation = true;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
elseif strcmpi(mesh,'H2H0')
    % LBB-Bedingung sollte verletzt werden (spurious oder checkerboard-Moden?)
    abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityVeryCoarse.inp');
    solidObject.elementDisplacementType = 'mixedSC';
    solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.flagHistoryFields = true;
    solidObject.mixedFEObject.condensation = false;
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.typeShapeFunctionData = 0;
elseif strcmpi(mesh,'H1H1')
    abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.flagHistoryFields = true;
    solidObject.mixedFEObject.continuousShapeFunctions = false;
    solidObject.mixedFEObject.condensation = true;
elseif strcmpi(mesh,'H2H2')
    % LBB-Bedingung sollte verletzt werden (spurious oder checkerboard-Moden?)
    abaqusMeshData = abaqusInputFileConverter('LShapeH2SerendipityVeryCoarse.inp');
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.flagHistoryFields = true;
    solidObject.mixedFEObject.condensation = true;
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.typeShapeFunctionData = 2;    
end
solidObject.meshObject.nodes = abaqusMeshData.qR;
if strcmpi(setupObject.integrator,'LinearImplicit')
    solidObject.meshObject.nodes = [solidObject.meshObject.nodes zeros(size(solidObject.meshObject.nodes))];
end
solidObject.meshObject.edof = abaqusMeshData.edof;
% solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.type = 'complex';

% solidObject.materialObject.name = 'ANN';
% solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
% solidObject.materialObject.name = 'MooneyRivlinVol2';
% solidObject.materialObject.name = 'MooneyRivlinModified';
solidObject.materialObject.name = 'MooneyRivlinModifiedVol2';
% solidObject.materialObject.name = 'SaintVenant';
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
solidObject.materialObject.lambda = 10000;
solidObject.materialObject.mu = mu;
if strcmpi(setupObject.integrator,'LinearImplicit')
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.rhoLinearImplicit = 100;
else
    solidObject.materialObject.rho = 100;
end
solidObject.dimension = 3;
% modified Mooney-Rivlin model (see Schroeder et al. 2011)
solidObject.materialObject.alpha = 42000;
solidObject.materialObject.beta = 84000;
solidObject.materialObject.gamma = 6*solidObject.materialObject.alpha + 12*solidObject.materialObject.beta;
solidObject.materialObject.epsilon1 = 100000;
solidObject.materialObject.epsilon2 = 10;

% solidObject.shapeFunctionObject.numberOfGausspoints = 27;

%% ANN
if strcmpi(solidObject.materialObject.name,'ANN')
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
end

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
nameFile = [solidObject.materialObject.name,solidObject.elementDisplacementType,solidObject.elementNameAdditionalSpecification,mesh,setupObject.integrator,'DT',num2str(setupObject.timeStepSize)];
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
try
    [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
    [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
end
strainEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
%externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = strainEnergy(1:length(kineticEnergy)) + kineticEnergy;% + externalEnergy(1:length(kineticEnergy));
if 1
    if strcmpi(setupObject.integrator,'DiscreteGradient')
        tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps)+1;
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
        % matlab2tikz(['totalAngularMomentum',nameFile,'.tikz'],'width','\figW','height','\figH')
    end
end
figure;
plot(timeVector, totalEnergy);
xlabel t; ylabel E; set(gca,'fontsize',20)
matlab2tikz(['energy',nameFile,'.tikz'],'width','\figW','height','\figH')
%
if 0
    figure;
    plot(timeVector(2:end), setupObject.newton.step','ro');
    xlabel t; ylabel Iterations; set(gca,'fontsize',20)
    % matlab2tikz(['iterations',nameFile,'.tikz'],'width','\figW','height','\figH')
end
%% postprocessing

myPath = 'simulationData/Plots'; % adjust the path

% plotTypes
plotTypes = {'totalEnergy', 'totalEnergyDiff', 'linearMomentum', 'totalLinearMomentum', 'angularMomentum', 'totalAngularMomentum', 'totalAngularMomentumDiff'};
 
% generate plots 
postprocessingPlots(dofObject,setupObject,bcTimeEnd,solidObject,mesh,plotTypes,myPath);
