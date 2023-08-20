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

% run('../../startUpMoofeKIT.m')
clear all;

neuronNameVector = [0,4,8,16,32,64,128];
computationTimeVector = zeros(1,numel(neuronNameVector));
coreCell = {'1 core', '24 cores'};
for coreCounter = 2%1:numel(coreCell)
    for neuronNumberCounter = 7%1:numel(neuronNameVector)

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'LShape';
        setupObject.saveObject.saveData = true;
        % setupObject.totalTimeSteps = 125;
        % setupObject.totalTimeSteps = 250;
        % setupObject.totalTimeSteps = 320;
        % setupObject.totalTime = 200;
        % setupObject.totalTimeSteps = 1;
        % setupObject.totalTime = 1;
        setupObject.totalTimeSteps = 250;
        setupObject.totalTime = 200;
        setupObject.plotObject.flag = false;
        % setupObject.plotObject.stress.component = -1;
        % setupObject.plotObject.view = [-0.4,90];
        setupObject.plotObject.view = [0,90];
        setupObject.plotObject.colorBarLimits = [0 400];
        setupObject.plotObject.border = 0.4;
        setupObject.plotObject.lineWidth = 2;
        
%         setupObject.newton.tolerance = 1e-5;
        setupObject.newton.tolerance = 5e-4;
        % setupObject.integrator = 'Endpoint';
        % setupObject.integrator = 'Midpoint';
        setupObject.integrator = 'DiscreteGradient';
        % 

        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
        % abaqus mesh
        mesh = 'H1';
%         mesh = 'H1H0'; 
%         mesh = 'H2H1';
        solidObject = solidClass(dofObject);
        if strcmpi(mesh,'H1')
        %     abaqusMeshData = abaqusInputFileConverter('LShapeH1VeryCoarse.inp');
        %     abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
            abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
            solidObject.shapeFunctionObject.order = 1;
            solidObject.shapeFunctionObject.numberOfGausspoints = 8;
            solidObject.elementDisplacementType = 'displacementSC';            
%             solidObject.mixedFEObject.typeShapeFunctionData = 0;
        elseif strcmpi(mesh,'H2H1')
            abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
            solidObject.shapeFunctionObject.order = 2;
            solidObject.shapeFunctionObject.numberOfGausspoints = 27;
            solidObject.mixedFEObject.typeShapeFunctionData = 1;
            solidObject.mixedFEObject.condensation = false;
            solidObject.elementDisplacementType = 'mixedSC';
            solidObject.mixedFEObject.condensation = true;        
        end
        poolObject = gcp('nocreate');
        if coreCounter == 1
            if ~isempty(poolObject)
                delete(poolObject);
            end
        elseif coreCounter == 2
            if isempty(poolObject)
                parpool(24)
            end
        end
        solidObject.meshObject.nodes = abaqusMeshData.qR;
        solidObject.meshObject.edof = abaqusMeshData.edof;
        % 
        solidObject.numericalTangentObject.computeNumericalTangent = false;
        solidObject.numericalTangentObject.showDifferences = false;
        if neuronNameVector(neuronNumberCounter) == 0
            %% GT
            solidObject.materialObject.name = 'MooneyRivlin';     
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
        else
            solidObject.materialObject.name = 'ANNMooneyRivlin';
            %% ANN
            % already trained ANN model
            solidObject.artificialNeuralNetworkObject.computeActivationFunction('softPlus');
            % solidObject.artificialNeuralNetworkObject.alpha = 5e2;
            neuronName = strcat('[',num2str(neuronNameVector(neuronNumberCounter)),']/');
            ANNlocationFolder = 'externalData/ANN/solid/normalized3/';
            solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = strcat(ANNlocationFolder,neuronName);
            solidObject.artificialNeuralNetworkObject.readData.biasNameCell = {'b1','b2'};
            solidObject.artificialNeuralNetworkObject.readData.weightNameCell = {'w1','w2'};
            solidObject.artificialNeuralNetworkObject.read('bias');
            solidObject.artificialNeuralNetworkObject.read('weight');
            solidObject.artificialNeuralNetworkObject.normalizationFactor = readmatrix(strcat(ANNlocationFolder,'normalizationFactor.txt'));
            solidObject.artificialNeuralNetworkObject.A = readmatrix(strcat(ANNlocationFolder,neuronName,'A.txt'));
        end
        solidObject.materialObject.rho = 100;
        solidObject.dimension = 3;

        FA = [256;512;768];
        bcTimeEnd = 5;
        neumannObject1 = neumannClass(dofObject);
        neumannObject1.loadType = 'deadLoad';
        neumannObject1.loadGeometry = 'area';
        neumannObject1.masterObject = solidObject;
        neumannObject1.loadVector = 1/9*FA;
        neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
        neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
        neumannObject1.projectionType = 'none';
        neumannObject1.timeFunction = @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
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
        neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;
        
        %% solver
        tic;
        %try
            dofObject = runNewton(setupObject,dofObject);
        %end
        computationTimeVector(coreCounter,neuronNumberCounter) = toc
    end
end

computationTimeVector
% computationTimeVector =  1.0e+03 *[   0.7673    1.4131    1.3150    1.4692    1.3712    1.4458    2.0202;
%                                         0.8370    0.8780    0.8910    0.9015    0.9009    0.9079        0.9444];

figure
neuronNameVector2 = 1:7;
bar(neuronNameVector2,computationTimeVector)
set(gca,'xticklabel',{'GT','4','8','16','32','64','128'});
legend('1 core','24 cores','location','north')
xlabel('GT or Number of Neurons')
ylabel('Computation time [s]')
set(gca,'fontsize',20)
nameFile = [solidObject.materialObject.name,setupObject.integrator,num2str(setupObject.timeStepSize)];
matlab2tikz(['computationTime',nameFile,'dofs',num2str(numel(solidObject.qN1)),'.tikz'],'width','\figW','height','\figH')