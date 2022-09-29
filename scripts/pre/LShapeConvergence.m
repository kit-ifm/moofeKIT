% LSHAPE Script for preprocessing a dynamic mechanical simulation.
% 
% FORMULATION
% Script generates a .mat file, which can be used to investigate the
% convergence order of the variables of a formulation.
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
% Felix Zaehringer

% run('../../startUpMoofeKIT.m');

computeInitialConfiguration = false;

if computeInitialConfiguration
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShapeConvergenceBasis';
setupObject.saveObject.saveData = true;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 5;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-4;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
solidObject = solidClass(dofObject);
solidObject.meshObject.nodes = abaqusMeshData.QNODE;
solidObject.meshObject.edof = abaqusMeshData.EDOF;
solidObject.materialObject.name = 'MooneyRivlin';
solidObject.elementDisplacementType = 'mixedSC';
%     solidObject.elementDisplacementType = 'displacement';
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
solidObject.mixedFEObject.orderShapeFunction = 1;

solidObject.materialObject.rho = 10;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;

bcTimeEnd = 5;
FA = [256;512;768];
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.typeOfLoad = 'deadLoad';
neumannObject1.masterObject = solidObject;
neumannObject1.forceVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.SUBSETS(3).EDOF;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.typeOfLoad = 'deadLoad';
neumannObject2.masterObject = solidObject;
neumannObject2.forceVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.SUBSETS(4).EDOF;

%% solver
warning off;
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
dofObject = runNewton(setupObject,dofObject);
end

clear;

variableTimeStepSizeVector = [0.025; 0.01; 0.005; 0.0025; 0.001; 0.0005];%[0.025;0.01;0.005;0.0025;0.001;0.0005];
for jj=1:size(variableTimeStepSizeVector,1)
    variableTimeStepSizeVector = [0.025; 0.01; 0.005; 0.0025; 0.001; 0.0005];
    % load Data of first steps
    load('LShapeConvergenceBasis', 'setupObject10', 'dofObject10');
    
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = strcat('LShapeConvergence', num2str(variableTimeStepSizeVector(jj)));
    setupObject.totalTimeSteps = 1/variableTimeStepSizeVector(jj);
    
    lastTimeStep = 1/variableTimeStepSizeVector(jj);
    setupObject.saveObject.saveData = true;
    setupObject.saveObject.saveForTimeSteps = lastTimeStep;
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    % setupObject.plotObject.stress.component = -1;
    setupObject.plotObject.view = [-0.4,90];
    setupObject.newton.tolerance = 1e-4;
    % setupObject.integrator = 'Endpoint';
%     setupObject.integrator = 'Midpoint';
    setupObject.integrator = 'DiscreteGradient';

    dofObject = dofClass;   % required object for dof and object handling

    %% continuum Objects
    % abaqus mesh
    abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
    solidObject = solidClass(dofObject);
    solidObject.meshObject.nodes = abaqusMeshData.QNODE;
    solidObject.meshObject.edof = abaqusMeshData.EDOF;
    solidObject.materialObject.name = 'MooneyRivlin';
%     solidObject.elementDisplacementType = 'displacementSC';
        solidObject.elementDisplacementType = 'mixedSC';
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
    solidObject.mixedFEObject.orderShapeFunction = 1;

    solidObject.materialObject.rho = 10;
    solidObject.dimension = 3;
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    
    
    solidObject.vN = dofObject10.listContinuumObjects{1}.vN1;
    solidObject.qN = dofObject10.listContinuumObjects{1}.qN1;
    solidObject.mixedFEObject.qN = dofObject10.listContinuumObjects{1}.mixedFEObject.qN1;
    
    %% solver
    warning off;
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    dofObject = runNewton(setupObject,dofObject);

    clear;
end



