flagNearlyIncompressible = false;
if flagNearlyIncompressible
    nameIncompressible = 'NearlyIncompressible';
else
    nameIncompressible = '';
end

%% setup (mandatory: setup and dofs)

% elementType = 'H1';
% elementType = 'H2';
% elementType = 'H1H1';
% elementType = 'H2H1';
% elementTypeCell = {'H1','H2H1'};
% elementTypeCell = {'H1H0','H2H1'};
elementTypeCell = {'H1','H1H0'};

% maxNumberOfElements = 2; % for H2H1
% maxNumberOfElements = 13; % for H2H1
maxNumberOfElements = 50;
convergenceVector = zeros(7,maxNumberOfElements);
convergenceVectorDOFs = zeros(7,maxNumberOfElements);
for model = 0% 0:6
    for elementTypeIndex = 1%1:numel(elementTypeCell)
        for nel = 6%1:maxNumberOfElements
            setupObject = setupClass;
            setupObject.saveObject.fileName = 'cooksMembrane';
            setupObject.saveObject.saveData = false;
            setupObject.totalTime = 1;
            setupObject.plotObject.flag = false;
            % setupObject.plotObject.stress.component = 11;
            setupObject.newton.tolerance = 1e-4;
            setupObject.integrator = 'Endpoint';

            dofObject = dofClass;   % required object for dof and object handling
            %% continuum Objects
            % [solidObject.nodes,solidObject.edof,edofNeumann] = linearTetrahedralCooksMembrane(nel,nel,nel,48,44,16,10);
            % FIXME: negative determinants for quadratic 27-node Element
            % [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = triquadraticCooksMembrane(nel,nel,nel,48,44,16,10);
            % FIXME: negative determinants for quadratic 10-node Element
            % [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = quadraticTetrahedralCooksMembrane(nel,nel,nel,48,44,16,10);

            % abaqusMeshData = abaqusInputFileConverter('cookSerendipity2x2x1.inp');
            % abaqusMeshData = abaqusInputFileConverter('cookSerendipity4x4x1.inp');
            %         abaqusMeshData = abaqusInputFileConverter('cooksMembrane.inp');

            % solidObject.meshObject.nodes = abaqusMeshData.qR;
            % solidObject.meshObject.edof = abaqusMeshData.edof;
            solidObject = solidClass(dofObject);
            if strcmpi(elementTypeCell{elementTypeIndex},'H1')
                solidObject.shapeFunctionObject.order = 1;
                solidObject.shapeFunctionObject.numberOfGausspoints = 8;
                solidObject.elementDisplacementType = 'displacementSC';
                setupObject.totalTimeSteps = 8;
            elseif strcmpi(elementTypeCell{elementTypeIndex},'H2')
                solidObject.shapeFunctionObject.order = 2;
                solidObject.shapeFunctionObject.numberOfGausspoints = 27;
                solidObject.elementDisplacementType = 'displacementSC';
                setupObject.totalTimeSteps = 4;
            elseif strcmpi(elementTypeCell{elementTypeIndex},'H1H1')
                solidObject.shapeFunctionObject.order = 1;
                solidObject.shapeFunctionObject.numberOfGausspoints = 8;
                solidObject.mixedFEObject.typeShapeFunctionData = 1;
                solidObject.elementDisplacementType = 'mixedSC';
                solidObject.mixedFEObject.condensation = true;
                setupObject.totalTimeSteps = 4;
            elseif strcmpi(elementTypeCell{elementTypeIndex},'H1H0')
                solidObject.shapeFunctionObject.order = 1;
                solidObject.shapeFunctionObject.numberOfGausspoints = 8;
                solidObject.mixedFEObject.typeShapeFunctionData = 0;
                solidObject.elementDisplacementType = 'mixedSC';
%                 solidObject.elementNameAdditionalSpecification = 'NonCascadeCGc';
%                 solidObject.elementNameAdditionalSpecification = 'InvariantCGc';
%                 solidObject.elementNameAdditionalSpecification = 'NonCascadeGc';
%                 solidObject.elementNameAdditionalSpecification = 'NonCascadeGJ';
                solidObject.elementNameAdditionalSpecification = 'NonCascadeIIIJ';

                % solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.condensation = false;
                solidObject.numericalTangentObject.computeNumericalTangent = true;
                solidObject.numericalTangentObject.type = 'complex';
                setupObject.totalTimeSteps = 8;
%                 setupObject.totalTimeSteps = 20;
            elseif strcmpi(elementTypeCell{elementTypeIndex},'H2H1')
                solidObject.shapeFunctionObject.order = 2;
                solidObject.shapeFunctionObject.numberOfGausspoints = 27;
                solidObject.mixedFEObject.typeShapeFunctionData = 1;
                % solidObject.elementDisplacementType = 'mixedSC';
                % solidObject.mixedFEObject.condensation = true;
                solidObject.elementDisplacementType = 'mixedSC';
                solidObject.elementNameAdditionalSpecification = 'InvariantCGc';
                solidObject.mixedFEObject.condensation = false;
                solidObject.numericalTangentObject.computeNumericalTangent = true;
                solidObject.numericalTangentObject.type = 'complex';
                setupObject.totalTimeSteps = 4;
            end
            solidObject.numericalTangentObject.computeNumericalTangent = true;
            if strcmpi(elementTypeCell{elementTypeIndex},'H1') || strcmpi(elementTypeCell{elementTypeIndex},'H1H1') || strcmpi(elementTypeCell{elementTypeIndex},'H1H0')
                [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 1, false);
%                 [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,1,48,44,16,4);
            elseif strcmpi(elementTypeCell{elementTypeIndex},'H2') || strcmpi(elementTypeCell{elementTypeIndex},'H2H1')
                [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 2, true);
            end
            %
            if model == 0
                solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
            else
                solidObject.materialObject.name = 'ANN';
                %% ANN
                % already trained ANN model
                solidObject.artificialNeuralNetworkObject.computeActivationFunction('softPlus');
                neuronNumberCounter = 2*2^model;
                % solidObject.artificialNeuralNetworkObject.alpha = 5e2;
                neuronName = strcat('[',num2str(neuronNumberCounter),']/');
                if ~flagNearlyIncompressible
                    ANNlocationFolder = 'externalData/ANN/solid/normalized3/';
                else
                    ANNlocationFolder = 'externalData/ANN/solid/normalized3Incompressible/';
                end
                solidObject.artificialNeuralNetworkObject.readData.locationFolderANN = strcat(ANNlocationFolder,neuronName);
                solidObject.artificialNeuralNetworkObject.readData.biasNameCell = {'b1','b2'};
                solidObject.artificialNeuralNetworkObject.readData.weightNameCell = {'w1','w2'};
                solidObject.artificialNeuralNetworkObject.read('bias');
                solidObject.artificialNeuralNetworkObject.read('weight');
                solidObject.artificialNeuralNetworkObject.normalizationFactor = readmatrix(strcat(ANNlocationFolder,'normalizationFactor.txt'));
                solidObject.artificialNeuralNetworkObject.A = readmatrix(strcat(ANNlocationFolder,neuronName,'A.txt'));
            end
            %
            if ~flagNearlyIncompressible
                bulk    = 5209;                                                 % K
                shear   = 997.5;                                                % G
                mu      = shear;                                                % second lame parameter
                a      = 5/6*mu;                                               % a
                b      = 1/6*mu;                                               % b
                c       = 10000;
            else
                a=126;
                b=252;
                c=81512;
            end
            solidObject.materialObject.a = a;
            solidObject.materialObject.b = b;
            solidObject.materialObject.c = c;
            solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
            solidObject.materialObject.rho = 0;
            solidObject.dimension = 3;

            dirichletObject = dirichletClass(dofObject);
            dirichletObject.nodeList = unique(bounEDOFs.SX1(:));
%             dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
            dirichletObject.nodalDof = 1:3;
            dirichletObject.masterObject = solidObject;
            dirichletObject.timeFunction = str2func('@(t) 0');

            neumannObject = neumannClass(dofObject);
            neumannObject.loadGeometry = 'area';
            neumannObject.loadType = 'deadLoad';
            neumannObject.masterObject = solidObject;
            if ~flagNearlyIncompressible
%                 neumannObject.loadVector = [0;500;0];
                neumannObject.loadVector = [0;0;200];
            else
                %             neumannObject.loadVector = [0;200;0];
                neumannObject.loadVector = [0;0;100];
%                 neumannObject.loadVector = [0;0;200];
            end
            neumannObject.projectionType = 'none';
            neumannObject.timeFunction = str2func('@(t) t');
%             neumannObject.meshObject.edof = edofNeumann;
            %         neumannObject.meshObject.edof = abaqusMeshData.subsets(4).edof;
            neumannObject.meshObject.edof = bounEDOFs.SX2;

            %% solver
            % try
                dofObject = runNewton(setupObject,dofObject);
                convergenceVector(model+1,nel) = max(solidObject.qN1(:,3))
                convergenceVectorDOFs(model+1,nel) = numel(solidObject.meshObject.nodes)
            % end
            % plot(solidObject,setupObject)
%             figure(99);
%             semilogy(1:size(setupObject.newton.normR(end-1,1:end),2),setupObject.newton.normR(end-1,1:end),'r-o','linewidth',2);
%             hold on
        end
    end
end
xlabel('Newton iteration'); ylabel('norm R'); set(gca,'fontsize',20)
legend('GT','ANN','Location','southwest')
matlab2tikz(['iterations',elementType,nameIncompressible,solidObject.elementDisplacementType,'.tikz'],'width','\figW','height','\figH')

convergenceVector = convergenceVector - max(solidObject.qR(:,3));
figure
% xValues = squeeze(convergenceVectorDOFs(1,:,:))';
% yValues = squeeze(convergenceVector(1,:,:))';
plot(convergenceVectorDOFs',convergenceVector','LineWidth',2)
% plot(convergenceVectorDOFs(2,2:end)',convergenceVector(2,2:end)','LineWidth',2)

xlabel('degrees of freedom')
ylabel('Displacement u_y^A [m]')
set(gca,'fontsize',20)
% legend('ANN-M-H2H1','Location','southeast')
% legend('GT-M-H2H1','Location','southeast')
legend('GT-H1','ANN-H1','Location','southeast')
matlab2tikz(['convergencePlot',elementTypeCell{elementTypeIndex},nameIncompressible,solidObject.materialObject.name,solidObject.elementDisplacementType,'.tikz'],'width','\figW','height','\figH')

if 0
        scrsz = get(groot,'screensize');
        bplot = scrsz(3)*0.2;
        hplot = scrsz(4)*0.4;
        xplot = scrsz(3)*[0.01,0.32,0.63];
        yplot = scrsz(4)*[0.5,0.1];
        figPlot = figure('outerposition',[xplot(3),yplot(2),bplot,hplot]);
        clf(figPlot)
        plot(solidObject,setupObject); 
        axis off; 
        axis equal;
        axis tight;
        view(66,28);
        if ~flagNearlyIncompressible
            caxis([50 550]);
        else
            caxis([50 350]);
        end
        colorbar off
        saveName = ['stressPlot',elementTypeCell{elementTypeIndex},nameIncompressible,solidObject.materialObject.name,solidObject.elementDisplacementType];
        print(figPlot,saveName,'-depsc')
        clf(figPlot)
        colormap('parula')
        colorbar
        axis off
        if ~flagNearlyIncompressible
            caxis([50 550]);
        else
            caxis([50 350]);
        end
        matlab2tikz(['Colorbar',saveName,'.tikz'],'width','\figW','height','\figH')
end