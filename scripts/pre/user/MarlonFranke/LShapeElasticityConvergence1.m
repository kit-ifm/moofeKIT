%% Calculation convergence
variableTimeStepSizeVector = [0.5; 0.25; 0.1; 0.05; 0.025; 0.01; 0.005; 0.0025; 0.001; 0.0005; 0.00025; 0.0001; 0.00001];
% variableTimeStepSizeVector = [0.25; 0.1; 0.025; 0.0001; 0.00001];
for jj=1:size(variableTimeStepSizeVector,1)
    % load Data of first steps
    % load('lShapeLoadPeriodMidpointLShapeH1.mat');
    load('lShapeLoadPeriodMidpointLShapeH1.mat');
    % load('rotatingXElectroThermoMechanicsBasis.mat', 'setupObject10', 'dofObject10');

    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.totalTime = 1;
    lastTimeStep = setupObject.totalTime/variableTimeStepSizeVector(jj);
    setupObject.totalTimeSteps = lastTimeStep;
    setupObject.saveObject.saveData = true;
    setupObject.saveObject.saveForTimeSteps = lastTimeStep;
    % if mod(lastTimeStep, 1000) == 0
    %     setupObject.saveObject.saveForTimeSteps = 1000:1000:lastTimeStep;
    % end
    setupObject.saveObject.saveForTimeSteps = lastTimeStep;
    setupObject.saveObject.exportLogFile = true;
    setupObject.plotObject.flag = false;
    setupObject.integrator = 'DiscreteGradient';
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.flagHistoryFields = true;
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.mixedFEObject.condensation = true;
    setupObject.newton.tolerance = 1e-4;
    setupObject.saveObject.fileName = strcat('lShape',solidObject.elementDisplacementType,solidObject.elementNameAdditionalSpecification,setupObject.integrator,'LShapeH1', num2str(variableTimeStepSizeVector(jj)));

    solidObject.mixedFEObject.qN = dofObject.listContinuumObjects{1}.mixedFEObject.qN1;
    solidObject.qN = dofObject.listContinuumObjects{1}.qN1;
    solidObject.vN = dofObject.listContinuumObjects{1}.vN1;    
    solidObject.mixedFEObject.qN = dofObject.listContinuumObjects{1}.mixedFEObject.qN1;
    
    clear neumannObject1 neumannObject2

    %% solver
    warning off;
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    try
        runNewton(setupObject,dofObject);
    end
end