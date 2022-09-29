% ROTATINGXELECTROTHERMOMECHANICS Script for preprocessing a dynamic
% coupled electro-thermo-mechanical simulation.
%
% REFERENCE
% In preparation
%
% SEE ALSO
% LShape
%
% CREATOR(S)
% Felix Zähringer, Marlon Franke

% run('../../startUpMoofeKIT.m');

computeInitialConfiguration = true;

armWidth = 1;
armLength = 3;
bodyThickness = 1;

if computeInitialConfiguration
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.saveData = false;
    setupObject.saveObject.fileName = 'rotatingXBasis';
    lastTimeStep = 120;
    setupObject.saveObject.saveForTimeSteps = 1:lastTimeStep;
    setupObject.totalTimeSteps = lastTimeStep;
    setupObject.totalTime = 0.6;
    setupObject.newton.tolerance = 1e-5;
    setupObject.newton.maximumResiduumNorm = 1e30;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.view = [-36.632241204539113,25.552647559456744];
    setupObject.integrator = 'DiscreteGradient';

    dofObject = dofClass;   % required object for dof and object handling

    %% continuum Objects
    solidObject = solidClass(dofObject);

    order = 2;
    serendipity = true;
    orderMixed = 1;
    numberOfGausspoints = 27;
    [nodes, edof] = meshRotatingX(armLength, armWidth, bodyThickness, 2, 1, 1, order, serendipity);
    solidObject.meshObject.nodes = nodes;
    solidObject.meshObject.edof = edof;

    solidObject.materialObject.name = 'MooneyRivlin';
    solidObject.elementDisplacementType = 'mixedSC';
%     solidObject.elementDisplacementType = 'displacementSC';
    solidObject.materialObject.rho = 1000;                                         % mass-density
    solidObject.materialObject.a = 25000;                                         % a
    solidObject.materialObject.b = 50000;                                         % b
    solidObject.materialObject.c = 500000;                                          % e
    solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);   % d
        
    solidObject.dimension = 3;
    solidObject.shapeFunctionObject.order = order;
    solidObject.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints;
    solidObject.mixedFEObject.condensation = true;
    solidObject.mixedFEObject.typeShapeFunctionData = orderMixed;

    % initial velocity
    vN = zeros(size(solidObject.meshObject.nodes));
    omega = [0, 0, 4];
    for ii=1:size(vN, 1)
        vN(ii,1:3) = cross(solidObject.meshObject.nodes(ii,1:3), omega);
    end
    solidObject.vN = vN;
    
    %% solver
    warning off;
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    runNewton(setupObject,dofObject);

    % clear;
end

% plot total energy
timeVector = getTime(dofObject.postDataObject, setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject, setupObject);
strainEnergy = getEnergy(dofObject.postDataObject, dofObject, setupObject, 'strainEnergy');
externalEnergy = getEnergy(dofObject.postDataObject, dofObject, setupObject, 'externalEnergy');
totalEnergy = kineticEnergy + strainEnergy;

figure;
plot(timeVector, totalEnergy);
xlim([timeVector(1), timeVector(end)]);
ylim([min(totalEnergy), max(totalEnergy)]);
xlabel('$t$ [s]');
ylabel('$E_{n+1}$ [J]');
figure;
totalEnergyN = totalEnergy(1:end-1);
totalEnergyDiff = totalEnergy(2:end) - totalEnergyN;
plot(timeVector(2:end), totalEnergyDiff);
xlabel('$t$ [s]');
ylabel('$E_{n+1}-E_{n}$ [J]');

%% Calculation convergence
computeConvergence = false;
if computeConvergence
    variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
    for jj=1:size(variableTimeStepSizeVector,1)
        variableTimeStepSizeVector = [0.001; 0.0005; 0.00005];
        % load Data of first steps
        load('rotatingXBasis.mat', 'setupObject10', 'dofObject10');
        
        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = strcat('rotatingXDT', num2str(variableTimeStepSizeVector(jj)));
        setupObject.totalTime = 0.1;
        lastTimeStep = setupObject.totalTime/variableTimeStepSizeVector(jj);
        
        setupObject.totalTimeSteps = lastTimeStep;
        
        setupObject.saveObject.saveData = true;
        setupObject.saveObject.saveForTimeSteps = lastTimeStep;
        if mod(lastTimeStep, 1000) == 0
            setupObject.saveObject.saveForTimeSteps = 1000:1000:lastTimeStep;
        end
        setupObject.saveObject.exportLogFile = true;
        setupObject.newton.tolerance = 5e-3;
        setupObject.newton.maximumResiduumNorm = 1e30;
%         setupObject.usePreconditioning = true;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.makeMovie = false;
        % setupObject.plotObject.postPlotType = 'stress';
        % setupObject.plotObject.postPlotType = 'none';
        % setupObject.plotObject.stress.component = -1;
        setupObject.plotObject.view = [-36.632241204539113,25.552647559456744];
        % setupObject.integrator = 'Endpoint';
        % setupObject.integrator = 'Midpoint';
        setupObject.integrator = 'DiscreteGradient';
        
%         if variableTimeStepSizeVector(jj) < 1e-4
%             setupObject.integrator = 'Midpoint';
%         else
%             setupObject.integrator = 'DiscreteGradient';
%         end
    
        dofObject = dofClass;   % required object for dof and object handling
    
        %% continuum Objects
        armWidth = 1;
        armLength = 3;
        bodyThickness = 1;


        solidObject = solidClass(dofObject);
        [nodes, edof] = meshRotatingX(armLength, armWidth, bodyThickness, 6, 2, 2, 2, 1);
        solidObject.meshObject.nodes = nodes;
        solidObject.meshObject.edof = edof;

        solidObject.materialObject.name = 'MooneyRivlin';
        solidObject.elementDisplacementType = 'mixedSC';
        solidObject.materialObject.rho = 1000;                                         % mass-density
        solidObject.materialObject.a = 25000;                                         % a
        solidObject.materialObject.b = 50000;                                         % b
        solidObject.materialObject.c = 500000;                                          % e
        solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);   % d
        
        solidObject.dimension = 3;
        solidObject.shapeFunctionObject.order = 2;
        solidObject.shapeFunctionObject.numberOfGausspoints = 27;
        solidObject.mixedFEObject.condensation = true;
        solidObject.mixedFEObject.orderShapeFunction = 1;
    
        solidObject.mixedFEObject.qN = dofObject10.listContinuumObjects{1}.mixedFEObject.qN1;
        solidObject.qN = dofObject10.listContinuumObjects{1}.qN1;
        solidObject.vN = dofObject10.listContinuumObjects{1}.vN1;

        %% solver
        warning off;
        parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
        parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
        runNewton(setupObject,dofObject);
        
        clear;
    end
end