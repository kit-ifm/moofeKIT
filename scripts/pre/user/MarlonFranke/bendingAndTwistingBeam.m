%
% FORMULATION
% Different formulations like standard displacement-based and mixed en-
% hanced assumed strain (eas) and different material models can be chosen.
%
% REFERENCE
% https://doi.org/10.1002/nme.5138
%
% SEE ALSO
% cooksMembrane,
% LShapeElectroThermo
%
% CREATOR(S)
% Marlon Franke

clear all;

% run('../../startUpMoofeKIT.m')

%% flag bending or twisting
flagBendingTwisting = 'bending';
% flagBendingTwisting = 'twisting';

flagFreeFlight = false;

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
% setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';
setupObject.integrator = 'LinearImplicit';
if strcmpi(flagBendingTwisting,'bending')
    setupObject.saveObject.fileName = 'bendingBeam';
    setupObject.plotObject.savePlot.name = 'bendingBeam';
    if strcmpi(setupObject.integrator,'LinearImplicit')
        setupObject.totalTimeSteps = 4000;
    else
        setupObject.totalTimeSteps = 40;
    end
    setupObject.totalTime = 2;
    % setupObject.plotObject.stress.component = -1;
    % setupObject.plotObject.view = [-0.4,90];
    setupObject.plotObject.view = [0,90];
    setupObject.plotObject.colorBarLimits = [0 4e5];
    setupObject.newton.tolerance = 5e-4;
    % setupObject.plotObject.xminxmaxset = [-5 5 -1 7 -4 4];
elseif strcmpi(flagBendingTwisting,'twisting')
    setupObject.saveObject.fileName = 'twistingBeam';
    setupObject.plotObject.savePlot.name = 'twistingBeam';
    if strcmpi(setupObject.integrator,'LinearImplicit')
        setupObject.totalTimeSteps = 400;
    else
        setupObject.totalTimeSteps = 40;
    end
    setupObject.totalTime = 0.4;
    % setupObject.plotObject.stress.component = -1;
    % setupObject.plotObject.view = [-0.4,90];
    setupObject.plotObject.view = [0,90];
    setupObject.plotObject.colorBarLimits = [0 2e6];
    setupObject.newton.tolerance = 5e-4;
    % setupObject.plotObject.xminxmaxset = [-3 3 -1 8 -4 4];
end
setupObject.saveObject.saveData = false;
setupObject.plotObject.flag = true;
setupObject.plotObject.savePlot.flag = false;
% setupObject.plotObject.steps = 100;
% setupObject.plotObject.keepFormerPlots = true;
setupObject.plotObject.savePlot.flagPlotColorbar = false;
setupObject.plotObject.border = 0.4;
setupObject.plotObject.lineWidth = 2;
setupObject.plotObject.colorscheme = 'turbo';
setupObject.plotObject.savePlot.type = '-dpng';
% setupObject.plotObject.savePlot.type = '-depsc';
% setupObject.newton.tolerance = 1e-6; % not good for ANN model<
% setupObject.newton.tolerance = 1e-8;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
mesh = 'H1H0';% mesh = 'H1';
% mesh = 'H2H1';% mesh = 'H2';
if strcmpi(setupObject.integrator,'LinearImplicit')
    solidObject = solidVelocityClass(dofObject);
else
    solidObject = solidClass(dofObject);
end
solidObject.dimension = 3;
if strcmpi(mesh,'H1') || strcmpi(mesh,'H1H0')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,6,1,2,12,2,1,false);
    % abaqusMeshData = abaqusInputFileConverter('LShapeH1Medium.inp');
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    % solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.mixedFEObject.typeShapeFunctionData = 0;
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.flagHistoryFields = true;
    % solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.mixedFEObject.condensation = true;
elseif strcmpi(mesh,'H2') || strcmpi(mesh,'H2H1')
    [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,6,1,2,12,2,2,true);
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    % solidObject.elementNameAdditionalSpecification = 'pHGJ';
    solidObject.flagHistoryFields = true;
    solidObject.shapeFunctionObject.order = 2;
    solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
    solidObject.mixedFEObject.condensation = false;
end
%% geometry adjustments
% solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1)+0.5;
solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2)+3;
% solidObject.meshObject.nodes(:,3) = solidObject.meshObject.nodes(:,1)+0.5;
solidObject.rotate(2,5.2*pi/180);

% solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.type = 'complex';
% solidObject.numericalTangentObject.showDifferences = true;

%% initial velocity
if strcmpi(flagBendingTwisting,'bending')
    if ~flagFreeFlight
        solidObject.vN = [5/3*solidObject.meshObject.nodes(:,2), zeros(size(solidObject.meshObject.nodes, 1),2)]; % Initialize velocity for each node
    else
        solidObject.vN = [5/3*solidObject.meshObject.nodes(:,2)+5/3*(solidObject.meshObject.nodes(:,2)-6), zeros(size(solidObject.meshObject.nodes, 1),2)]; % Initialize velocity for each node
    end
elseif strcmpi(flagBendingTwisting,'twisting')
    solidObject.vN = 100*sin(pi*solidObject.meshObject.nodes(:,2)/12).*[solidObject.meshObject.nodes(:,3), zeros(size(solidObject.meshObject.nodes, 1),1), -solidObject.meshObject.nodes(:,1)]; % Initialize velocity for each node
else
    error('flagBendingTwisting need to be either bending or twisting')
end
if strcmpi(setupObject.integrator,'LinearImplicit')
    if strcmpi(flagBendingTwisting,'bending')
        solidObject.meshObject.nodes = [solidObject.meshObject.nodes [5/3*solidObject.meshObject.nodes(:,2), zeros(size(solidObject.meshObject.nodes, 1),2)]];
    elseif strcmpi(flagBendingTwisting,'twisting')
        solidObject.meshObject.nodes = [solidObject.meshObject.nodes 100*sin(pi*solidObject.meshObject.nodes(:,2)/12).*[solidObject.meshObject.nodes(:,3), zeros(size(solidObject.meshObject.nodes, 1),1), -solidObject.meshObject.nodes(:,1)]];
    end
    solidObject.vN = [solidObject.vN, 0*solidObject.vN];
end

%
% solidObject.materialObject.name = 'MooneyRivlinModified';
% solidObject.materialObject.name = 'MooneyRivlinVol2';
solidObject.materialObject.name = 'MooneyRivlinModifiedVol2';
% solidObject.materialObject.name = 'MooneyRivlin';
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

% modified Mooney-Rivlin model (see Schroeder et al. 2011)
solidObject.materialObject.alpha = 42000;
solidObject.materialObject.beta = 84000;
solidObject.materialObject.gamma = 6*solidObject.materialObject.alpha + 12*solidObject.materialObject.beta;
solidObject.materialObject.epsilon1 = 100000;
solidObject.materialObject.epsilon2 = 10;
solidObject.materialObject.c = 1e5;
% solidObject.materialObject.c = 1e6; % kommt das gleiche raus
solidObject.materialObject.d = 6*solidObject.materialObject.alpha + 12*solidObject.materialObject.beta;

if strcmpi(setupObject.integrator,'LinearImplicit')
    solidObject.materialObject.rho = 0;
    solidObject.materialObject.rhoLinearImplicit = 225;
    % solidObject.materialObject.rhoLinearImplicit = 1100;
else
    solidObject.materialObject.rho = 225;
    % solidObject.materialObject.rho = 1100;
end
solidObject.dimension = 3;

if ~flagFreeFlight
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,2)==0);
    if strcmpi(setupObject.integrator,'LinearImplicit')
        dirichletObject.nodalDof = 1:6;
    else
        dirichletObject.nodalDof = 1:3;
    end
    dirichletObject.masterObject = solidObject;
    dirichletObject.timeFunction = str2func('@(t) 0');
end

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
bcTimeEnd = 0;
nameFile = [solidObject.materialObject.name,solidObject.elementDisplacementType,solidObject.elementNameAdditionalSpecification,mesh,setupObject.integrator,'DT',num2str(setupObject.timeStepSize)];
timeVector = getTime(dofObject.postDataObject,setupObject);
if strcmpi(setupObject.integrator,'LinearImplicit')
    kineticEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'kineticEnergy');
else
    kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
end
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
if strcmpi(flagBendingTwisting,'bending')
    myPath = 'simulationData/PlotsBendingBeam'; % adjust the path
elseif strcmpi(flagBendingTwisting,'twisting')
    myPath = 'simulationData/PlotsTwistingBeam'; % adjust the path
else
    error('flagBendingTwisting need to be either bending or twisting')
end

% plotTypes
plotTypes = {'totalEnergy', 'totalEnergyDiff', 'linearMomentum', 'totalLinearMomentum', 'angularMomentum', 'totalAngularMomentum', 'totalAngularMomentumDiff'};
 
% generate plots 
postprocessingPlots(dofObject,setupObject,bcTimeEnd,solidObject,mesh,plotTypes,myPath);
