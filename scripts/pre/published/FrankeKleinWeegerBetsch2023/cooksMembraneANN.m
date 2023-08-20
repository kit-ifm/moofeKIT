clear all;

flagNearlyIncompressible = true;
if flagNearlyIncompressible
    nameIncompressible = 'NearlyIncompressible';
else
    nameIncompressible = '';
end

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembrane';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
% setupObject.plotObject.stress.component = 11;
setupObject.newton.tolerance = 1e-6;
setupObject.integrator = 'Endpoint';

dofObject = dofClass;   % required object for dof and object handling
%% continuum Objects
solidObject = solidClass(dofObject);
% mesh = 'H1'; % 10 load steps needed
% mesh = 'H2';
mesh = 'H1H0';
% mesh = 'H1H1';
% mesh = 'H2H1';
% mesh = 'H2H2';
nel = 4;
if strcmpi(mesh,'H1')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 1, false);
%     [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,1,48,44,16,4);
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
elseif strcmpi(mesh,'H1H0')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 1, false);
%     [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,1,48,44,16,4);
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 0;
elseif strcmpi(mesh,'H1H1')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 1, false);
%     [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,1,48,44,16,4);
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
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
elseif strcmpi(mesh,'H2H2')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof,bounEDOFs] = meshGeneratorCooksMembrane(nel, 1, nel, 2, true);
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.typeShapeFunctionData = 2;
end
solidObject.mixedFEObject.condensation = false;
solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.epsilon = 1e-6;
% solidObject.numericalTangentObject.type = 'centralDifferences';
solidObject.numericalTangentObject.type = 'complex';
solidObject.numericalTangentObject.showDifferences = false;
%
% solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
solidObject.materialObject.name = 'ANN';
% solidObject.elementDisplacementType = 'displacementSC';
solidObject.elementDisplacementType = 'mixedSC';
% solidObject.elementNameAdditionalSpecification = 'NonCascadeCGc';
% solidObject.elementNameAdditionalSpecification = 'InvariantCGc';
% solidObject.elementNameAdditionalSpecification = 'NonCascadeGc';
% solidObject.elementNameAdditionalSpecification = 'NonCascadeGJ';
solidObject.elementNameAdditionalSpecification = 'NonCascadeIIIJ';
if ~flagNearlyIncompressible
    bulk    = 5209;                                                 % K
    shear   = 997.5;                                                % G
    mu      = shear;                                                % second lame parameter
    a      = 5/6*mu;                                               % a
    b      = 1/6*mu;                                               % b
    c       = 10000;
else
    a = 126;
    b = 252;
    c = 81512;
end
solidObject.materialObject.a = a;
solidObject.materialObject.b = b;
solidObject.materialObject.c = c;
solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
solidObject.materialObject.rho = 0;
solidObject.dimension = 3;

%% ANN
% already trained ANN model
solidObject.artificialNeuralNetworkObject.computeActivationFunction('softPlus');
neuronNumberCounter = 8;
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

dirichletObject = dirichletClass(dofObject);
% dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
dirichletObject.nodeList = unique(bounEDOFs.SX1(:));
dirichletObject.nodalDof = 1:3;
dirichletObject.masterObject = solidObject;
dirichletObject.timeFunction = str2func('@(t) 0');

neumannObject = neumannClass(dofObject);
neumannObject.loadType = 'deadLoad';
neumannObject.masterObject = solidObject;
if ~flagNearlyIncompressible
    neumannObject.loadVector = [0;0;500];
%     neumannObject.loadVector = [0;0;100];
else
%     neumannObject.loadVector = [0;0;100];
    neumannObject.loadVector = [0;0;100];
end
neumannObject.loadGeometry = 'area';
% neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
% neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject.projectionType = 'none';
neumannObject.timeFunction = str2func('@(t) t');
% neumannObject.meshObject.edof = edofNeumann;
% neumannObject.meshObject.edof = abaqusMeshData.subsets(4).edof;
neumannObject.meshObject.edof = bounEDOFs.SX2;

%% solver
dofObject = runNewton(setupObject,dofObject);
