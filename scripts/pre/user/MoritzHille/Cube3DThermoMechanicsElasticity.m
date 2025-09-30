% 3D Cube Script for preprocessing a dynamic mechanical simulation.
%
% FORMULATION
% Different formulations like standard displacement-based and mixed en-
% hanced assumed strain (eas) and different material models can be chosen.
%
% REFERENCE
% https://doi.org/10.1007/BF00913408
%
%
% CREATOR(S)
% 08.02.2025 Tim Pleßke, Moritz Hille
%% color scheme for plots
color_scheme = {'#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00','#CC79A7'};

%% setup (mandatory: setup and dofs)
setupObject                     = setupClass;
setupObject.saveObject.fileName = 'Cube3DThermoElasticity';
setupObject.saveObject.saveData = false;
setupObject.newton.tolerance    = 1e-9;
setupObject.newton.maximumSteps = 500;

testmode = 'stress';
% testmode = 'temp';


% setupObject.integrator = 'Endpoint';
setupObject.integrator   = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

setupObject.plotObject.flag             = true;
setupObject.plotObject.postPlotType     = testmode;

if strcmp(testmode,'temp')
    setupObject.plotObject.colorBarLimits = [270 310];
    setupObject.totalTimeSteps = 700;
    setupObject.totalTime = 4;
else
    setupObject.plotObject.colorBarLimits = [0 600];
    setupObject.totalTimeSteps = 250;
    setupObject.totalTime = 6;
end

% setupObject.plotObject.stress.component = -1;
% setupObject.plotObject.view             = [20,80,80];
% setupObject.plotObject.view             = [10,10,10];
setupObject.plotObject.view               = [-20,45];
setupObject.plotObject.border             = 0.4;
setupObject.plotObject.lineWidth          = 2;
setupObject.plotObject.savePlot.flag      = true;
setupObject.plotObject.savePlot.name      = 'CubeGT';
setupObject.plotObject.savePlot.type      = '-depsc';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects

% mesh generator 3D cube
nelXYZ = 4;
solidThermoObject = solidThermoClass(dofObject);
lengthX           = 1; % length of brick in x-direction
lengthY           = 1; % length of brick in y-direction
lengthZ           = 1; % length of brick in z-direction
nelX              = nelXYZ; % number of elements in x-direction
nelY              = nelXYZ; % number of elements in y-direction
nelZ              = nelXYZ; % number of elements in z-direction
serendipity       = false;

solidThermoObject.dimension                               = 3;
solidThermoObject.shapeFunctionObject.order               = 1;
solidThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
solidThermoObject.mixedFEObject.typeShapeFunctionData     = 1;

[solidThermoObject.meshObject.nodes, solidThermoObject.meshObject.edof, boundEDOF] = meshGeneratorCube(lengthX, lengthY, lengthZ, nelX, nelY, nelZ, solidThermoObject.shapeFunctionObject.order , serendipity);

% initial thermal field
nnodes=size(solidThermoObject.meshObject.nodes(:,1),1);
initialThermalField = 293.15*ones(nnodes,1);

% inital thermal field at the left side of the cube (testmode temp)
if strcmp(testmode,'temp')
    % initialThermalField([boundEDOF.SY1],1) = 313.15;
    % initialThermalField([boundEDOF.SY2],1) = 313.15;
    bcTimeEnd = 2;
    neumannObject1 = neumannClass(dofObject);
    neumannObject1.masterObject = solidThermoObject;
    neumannObject1.loadPhysics = 'thermal';
    neumannObject1.loadGeometry = 'area';
    neumannObject1.loadVector = 3;
    neumannObject1.projectionType = 'none';
    neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
    neumannObject1.meshObject.edof = boundEDOF.SY1;
    neumannObject2 = neumannClass(dofObject);
    neumannObject2.masterObject = solidThermoObject;
    neumannObject2.loadPhysics = 'thermal';
    neumannObject2.loadGeometry = 'area';
    neumannObject2.loadVector = -3;
    neumannObject2.projectionType = 'none';
    neumannObject2.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
    neumannObject2.meshObject.edof = boundEDOF.SY2;
end

solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, initialThermalField];

% material model and element type

% solidViscoObject.materialObject.name  = 'Hooke';
% solidThermoObject.materialObject.name = 'NeoHookeGENERIC';
% solidThermoObject.materialObject.name = 'MooneyRivlin';
% solidThermoObject.materialObject.name = 'MooneyRivlinWedge';
solidThermoObject.materialObject.name = 'MooneyRivlinGENERIC';
% solidThermoObject.materialObject.name = 'SaintVenant';
% solidThermoObject.materialObject.name = 'SaintVenantNumericalTangent';
% solidThermoObject.materialObject.name = 'Hyperelastic'; % 'Hooke';

% solidThermoObject.elementDisplacementType = 'displacement';
%   solidThermoObject.elementDisplacementType = 'dispL2';
solidThermoObject.elementDisplacementType = 'mixedCGJL2';
% solidThermoObject.elementDisplacementType = 'displacementSC';
% solidThermoObject.elementDisplacementType = 'eas';
% solidThermoObject.mixedFEObject.condensation = true;

% material data
E                                       = 2.26115e+03;
nu                                      = 3.95469e-01;
% solidThermoObject.materialObject.rho    = 10;                                                                               % mass density
solidThermoObject.materialObject.rho    = 0;                                                                               % mass density
solidThermoObject.materialObject.lambda = 5209;
solidThermoObject.materialObject.mu     = 997.5;
solidThermoObject.materialObject.a      = solidThermoObject.materialObject.mu/4;                                            % material parameter a
solidThermoObject.materialObject.b      = solidThermoObject.materialObject.mu/4;                                            % material parameter b
solidThermoObject.materialObject.c1     = solidThermoObject.materialObject.lambda - solidThermoObject.materialObject.mu;    % material parameter c1                                                                 % material parameter c1
solidThermoObject.materialObject.d1     = 2*(solidThermoObject.materialObject.a + 2*solidThermoObject.materialObject.b);    % material parameter d1
solidThermoObject.materialObject.c2     = solidThermoObject.materialObject.lambda;                                          % material parameter c2
solidThermoObject.materialObject.d2     = 2/3 * solidThermoObject.materialObject.mu;                                        % material parameter d2

solidThermoObject.materialObject.kappa  = 100;                                                                              % heat capacity
solidThermoObject.materialObject.beta   = 2.233*10^(-4);                                                                    % coupling parameter
solidThermoObject.materialObject.thetaR = 293.15;                                                                           % reference Temperature
solidThermoObject.materialObject.k0     = 10;                                                                               % thermal conductivity

% numerical tangent
solidThermoObject.numericalTangentObject.computeNumericalTangent = false;
solidThermoObject.numericalTangentObject.showDifferences = false;
solidThermoObject.numericalTangentObject.type = 'complex';

% boundary conditions (testmode stress)
bcTimeEnd = 2;
FA = [0;0;1];

if strcmp(testmode,'stress')
    solidThermoObject.materialObject.beta = 0;

    neumannObject1 = neumannClass(dofObject);
    neumannObject1.masterObject = solidThermoObject;
    neumannObject1.loadGeometry = 'area';
    neumannObject1.loadVector = -1000*FA;
    neumannObject1.projectionType = 'none';
    neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
    neumannObject1.meshObject.edof = boundEDOF.SZ2;
end

% dirichlet boundaries
dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.masterObject = solidThermoObject;
dirichletObject1.nodeList = 1;
dirichletObject1.nodalDof = 1;
dirichletObject1.timeFunction = str2func(['@(t, X) ', num2str(-lengthX/2)]);

dirichletObject3 = dirichletClass(dofObject);
dirichletObject3.masterObject = solidThermoObject;
dirichletObject3.nodeList = 1;
dirichletObject3.nodalDof = 2;
dirichletObject3.timeFunction = str2func(['@(t, Y) ', num2str(-lengthY/2)]);

dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.masterObject = solidThermoObject;
dirichletObject2.nodeList = boundEDOF.SZ1;
dirichletObject2.nodalDof = 3;
dirichletObject2.timeFunction = str2func(['@(t, Z) ', num2str(-lengthZ/2)]);


%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy

% generate required data
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
[linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = internalEnergy + kineticEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+3:end) - totalEnergy(tStartDiff+2:end-1);
totalEntropy = getElementData(dofObject.postDataObject,dofObject,setupObject,'entropy');
totalEntropyDiff = getElementData(dofObject.postDataObject,dofObject,setupObject,'deltaS');
figure;

% initialize folder to save plots

% base folder for saving plots
% baseFolder = '/Users/timplesske/Documents/MATLAB/Plots';
baseFolder = '/Users/Moritz Hille/Desktop/Nextcloud/MarkplatzMH/Forschung/GENERIC/Programmierung/Ergebnisse Matlab'; % Adjust the path
% baseFolder = '/Users/morit/Desktop/Ergebnisse Matlab'; % Adjust the path
% create folder name for saving plots
if solidThermoObject.numericalTangentObject.computeNumericalTangent
    folderName = sprintf('%s_%s_%s_%s_numTang', ...
        setupObject.saveObject.fileName, ...
        solidThermoObject.elementDisplacementType, ...
        solidThermoObject.materialObject.name, ...
        setupObject.integrator);
else
    folderName = sprintf('%s_%s_%s_%s', ...
        setupObject.saveObject.fileName, ...
        solidThermoObject.elementDisplacementType, ...
        solidThermoObject.materialObject.name, ...
        setupObject.integrator);
end

folder = fullfile(baseFolder, folderName);

% folder version
version = 0;
while exist(folder, 'dir')
    version = version + 1;
    folder = fullfile(baseFolder, sprintf('%s_V%d', folderName, version));
end

% create the main folder
mkdir(folder);

% define plot types and their descriptions
plotTypes = {'totalEnergy', 'totalEnergyDiff', 'linearMomentum', ...
    'totalLinearMomentum', 'angularMomentum', 'totalAngularMomentum', 'totalEntropy', 'totalEntropyDiff'};

% loop over all plot types
for i = 1:length(plotTypes)
    % create filenames for .fig and .tikz files
    if version > 0
        % for version > 1 , the version is added to the filename
        if solidThermoObject.numericalTangentObject.computeNumericalTangent
            figName = sprintf('%s_%s_%s_%s_numTang_%s_V%d.fig', ...
                setupObject.saveObject.fileName, ...
                solidThermoObject.elementDisplacementType, ...
                solidThermoObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i}, version);
        else
            figName = sprintf('%s_%s_%s_%s_%s_V%d.fig', ...
                setupObject.saveObject.fileName, ...
                solidThermoObject.elementDisplacementType, ...
                solidThermoObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i}, version);
        end
    else
        % for the first save (version 0), no version is added
        if solidThermoObject.numericalTangentObject.computeNumericalTangent
            figName = sprintf('%s_%s_%s_%s_numTang_%s.fig', ...
                setupObject.saveObject.fileName, ...
                solidThermoObject.elementDisplacementType, ...
                solidThermoObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i});
        else
            figName = sprintf('%s_%s_%s_%s_%s.fig', ...
                setupObject.saveObject.fileName, ...
                solidThermoObject.elementDisplacementType, ...
                solidThermoObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i});
        end
    end

    % full path for the .fig file
    figPath = fullfile(folder, figName);

    % full path for the .tikz file
    tikzPath = strrep(figPath, '.fig', '.tikz');

    % create the plot
    switch plotTypes{i}
        case 'totalEnergy'
            plot(timeVector, internalEnergy, 'color', color_scheme{1}, 'LineWidth', 0.9);
            hold on;
            plot(timeVector, kineticEnergy, 'color', color_scheme{5}, 'LineWidth', 0.9);
            plot(timeVector, totalEnergy, 'color', color_scheme{7}, 'LineWidth', 0.9);
            title(['internal, kinetic and total energy with $\beta = $', num2str(solidThermoObject.materialObject.beta)], 'Interpreter', 'latex');
            legend('internal energy', 'kinetic energy', 'total energy');

        case 'totalEnergyDiff'
            plot(timeVector(tStartDiff + 3:end), totalEnergyDiff, 'color', color_scheme{3}, 'LineWidth', 0.9);
            legend('total energy diff');
            title(['total energy diff with $\beta = $', num2str(solidThermoObject.materialObject.beta)], 'Interpreter', 'latex');

        case 'linearMomentum'
            plot(timeVector, linearMomentum);
            legend('linear momentum');
            title('linear momentum');

        case 'totalLinearMomentum'
            plot(timeVector, totalLinearMomentum);
            legend('total linear momentum');
            title('total linear momentum');

        case 'angularMomentum'
            plot(timeVector, angularMomentum);
            title('angular momentum');

        case 'totalAngularMomentum'
            plot(timeVector, totalAngularMomentum, 'color', color_scheme{1}, 'LineWidth', 0.9);
            title(['total angular momentum with $\beta = $', num2str(solidThermoObject.materialObject.beta)], 'Interpreter', 'latex');

        case 'totalEntropy'
            plot(timeVector, totalEntropy, 'color', color_scheme{1}, 'LineWidth', 0.9);
            title('total entropy');
            legend('total entropy')

        case 'totalEntropyDiff'
            plot(timeVector(2:end), totalEntropyDiff(2:end), 'color', color_scheme{3}, 'LineWidth', 0.9);
            legend('total entropy diff');
            title(['total entropy diff with $\beta = $', num2str(solidThermoObject.materialObject.beta)], 'Interpreter', 'latex');
    end

    % save the plot
    savefig(figPath);
    matlab2tikz('filename', tikzPath);
    if i < length(plotTypes)
        figure; % create new figure for next plot
    end
end