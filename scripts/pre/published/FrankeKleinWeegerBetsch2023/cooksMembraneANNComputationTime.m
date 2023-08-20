clear all;

%% continuum Objects
mesh = 'H1'; % 10 load steps needed
% mesh = 'H2';
% mesh = 'H2H1';
% nel = 96;
% nel = 48;
nel = 24;
neuronNameVector = [0,4,8,16,32,64,128];
computationTimeVector = zeros(1,numel(neuronNameVector));
coreCell = {'1 core', '24 cores'};
for coreCounter = 1:numel(coreCell)
    for neuronNumberCounter = 1:numel(neuronNameVector)
        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'cooksMembrane';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 4;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        % setupObject.plotObject.stress.component = 11;
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';
        dofObject = dofClass;   % required object for dof and object handling
        solidObject = solidClass(dofObject);
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
        if strcmpi(mesh,'H1')
%             [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 8, nel, 1, false);
            [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 4, nel, 1, false);
%             [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 2, nel, 1, false);
%             [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 1, false);
            %     [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,1,48,44,16,4);
            solidObject.shapeFunctionObject.order = 1;
            solidObject.shapeFunctionObject.numberOfGausspoints = 8;
        elseif strcmpi(mesh,'H2H1')
            [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 2, true);
            %     solidObject.meshObject.nodes = solidObject.meshObject.nodes +  [20    5   25];
            %         [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = triquadraticCooksMembrane(nel,nel,nel,48,44,16,10);
            % FIXME: negative determinants for quadratic 10-node Element
            %     [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = quadraticTetrahedralCooksMembrane(nel,nel,nel,48,44,16,10);
            % abaqusMeshData = abaqusInputFileConverter('cookSerendipity2x2x1.inp');
            % abaqusMeshData = abaqusInputFileConverter('cookSerendipity4x4x1.inp');
            %     abaqusMeshData = abaqusInputFileConverter('cooksMembrane.inp');
            %     solidObject.meshObject.nodes = abaqusMeshData.qR+[53,20,0];
            %     solidObject.meshObject.nodes = solidObject.meshObject.nodes +[0,0,10];
            %     solidObject.meshObject.edof = abaqusMeshData.edof;
            solidObject.shapeFunctionObject.order = 2;
            solidObject.shapeFunctionObject.numberOfGausspoints = 27;
            solidObject.mixedFEObject.typeShapeFunctionData = 1;
        end
        solidObject.mixedFEObject.condensation = true;

        solidObject.numericalTangentObject.computeNumericalTangent = false;
        solidObject.numericalTangentObject.showDifferences = false;
        %
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
        solidObject.materialObject.rho = 0;
        solidObject.dimension = 3;
        solidObject.elementDisplacementType = 'displacementSC';

        dirichletObject = dirichletClass(dofObject);
        % dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
        dirichletObject.nodeList = unique(bounEDOFs.SX1(:));
        dirichletObject.nodalDof = 1:3;
        dirichletObject.masterObject = solidObject;
        dirichletObject.timeFunction = str2func('@(t) 0');

        neumannObject = neumannClass(dofObject);
        neumannObject.loadType = 'deadLoad';
        neumannObject.loadGeometry = 'area';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVector = [0;0;500];
        neumannObject.projectionType = 'none';
        neumannObject.timeFunction = str2func('@(t) t');
        % neumannObject.meshObject.edof = edofNeumann;
        % neumannObject.meshObject.edof = abaqusMeshData.subsets(4).edof;
        neumannObject.meshObject.edof = bounEDOFs.SX2;

        %% solver
        tic;
        %try
            dofObject = runNewton(setupObject,dofObject);
        %end
        computationTimeVector(coreCounter,neuronNumberCounter) = toc
    end
end

if 1
    figure
    neuronNameVector2 = 1:7;
    bar(neuronNameVector2,computationTimeVector)
    set(gca,'xticklabel',{'GT','4','8','16','32','64','128'});
    legend('1 core','24 cores','location','north')
    xlabel('GT or Number of Neurons')
    ylabel('Computation time [s]')
    set(gca,'fontsize',20)
    nameFile = [solidObject.materialObject.name,setupObject.integrator,num2str(setupObject.timeStepSize)];
    matlab2tikz(['computationTime24Cores',nameFile,num2str(numel(solidObject.qN1)),'Elements','.tikz'],'width','\figW','height','\figH')
end

if 0
    loadsteps = 1:setupObject.totalTimeSteps;
    iterationNumberGT = [ 7 7 6 6];
    iterationNumberANN8 = [ 6 7 6 6];
    figure
    plot(loadsteps,iterationNumberGT,'ro',loadsteps,iterationNumberANN8,'bx')
    legend('GT','ANN','location','north')
    xlabel('load step')
    ylabel('Iterations')
    set(gca,'fontsize',20)
%     set(gca,'yticklabel',{'6','7'});
    nameFile = [solidObject.materialObject.name,setupObject.integrator,num2str(setupObject.timeStepSize)];
    matlab2tikz(['iterations',nameFile,num2str(numel(solidObject.qN1)),'Elements','.tikz'],'width','\figW','height','\figH')
end
