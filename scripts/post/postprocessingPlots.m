function [] = postprocessingPlots(dofObject,setupObject,bcTimeEnd,continuumObject,mesh,plotTypes,myPath)
%% postprocessing function
%
% Generates and saves the following plots:
% - Energy: totalEnergy, totalEnergyDiff
% - Momentum: linearMomentum, totalLinearMomentum, angularMomentum, totalAngularMomentum
% - Entropy: totalEntropy, totalEntropyDiff
%
% File and folder structure:
% - Creates a folder in 'myPath' named:
%   'setupObject.saveObject.fileName_solidObject.elementDisplacementType_ ...
%   solidObject.materialObject.name_setupObject.integrator_V{i}'
% - Saves '.fig' and '.tikz' files in 'myPath'.
% - Files follow the naming pattern '[folder_name]_[plotType]'.
% - If numerical tangents are computed, the following naming pattern is applied:
%   '[folder_name]_numTang_[plotType]'.
%
% Set your height, width, fontsize, and linewidth for your tikzpicture
% - For example:
%   \setlength{\figH}{5cm}  (command in ifmclasses)
%   \setlength{\figW}{7.5cm}  (command in ifmclasses)
%   \newcommand{\figFontSize}{12pt}
%   \newcommand{\figLineWidth}{1.2pt}
%
% 25.02.2025 Tim Pleßke

% color scheme for plots
color_scheme = {'#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00','#CC79A7'};

% generate required data
timeVector = getTime(dofObject.postDataObject,setupObject);
timeStepSize = setupObject.totalTime / setupObject.totalTimeSteps;


kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
if isa(continuumObject, 'solidClass')
    internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
else
    internalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
end
externalEnergy = getElementData(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');

totalEnergy = internalEnergy + kineticEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+3:end) - totalEnergy(tStartDiff+2:end-1);

[linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);

totalAngularMomentumDiff = totalAngularMomentum(tStartDiff+3:end) - totalAngularMomentum(tStartDiff+2:end-1);

totalEntropy = getElementData(dofObject.postDataObject,dofObject,setupObject,'entropy');
totalEntropyDiff = getElementData(dofObject.postDataObject,dofObject,setupObject,'deltaS');
figure;

% initialize folder to save plots

% base folder for saving plots
baseFolder = myPath;

% create folder name for saving plots
if continuumObject.numericalTangentObject.computeNumericalTangent
    tangentName = '_numTang';
else
    tangentName = '';
end
if ~isempty(continuumObject.elementNameAdditionalSpecification)
    elementNameAdditionalSpecification = strcat('_',continuumObject.elementNameAdditionalSpecification);
else
    elementNameAdditionalSpecification = '_';
end
folderName = sprintf('%s_%s_%s%s_%s_%s', ...
    setupObject.saveObject.fileName, ...
    continuumObject.materialObject.name, ...
    continuumObject.elementDisplacementType, ...
    elementNameAdditionalSpecification, ...
    mesh, ...
    setupObject.integrator,...
    tangentName);

folder = fullfile(baseFolder, folderName);

% folder version
version = 0;
while exist(folder, 'dir')
    version = version + 1;
    folder = fullfile(baseFolder, sprintf('%s_V%d', folderName, version));
end

% create the main folder
mkdir(folder);

% loop over all plot types
for i = 1:length(plotTypes)
    % create filenames for .fig and .tikz files
    if version > 0
        % for version > 1 , the version is added to the filename
        if continuumObject.numericalTangentObject.computeNumericalTangent
            figName = sprintf('%s_%s_%s_%s_numTang_%s_V%d.fig', ...
                setupObject.saveObject.fileName, ...
                continuumObject.elementDisplacementType, ...
                continuumObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i}, version);
        else
            figName = sprintf('%s_%s_%s_%s_%s_V%d.fig', ...
                setupObject.saveObject.fileName, ...
                continuumObject.elementDisplacementType, ...
                continuumObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i}, version);
        end
    else
        % for the first save (version 0), no version is added
        if continuumObject.numericalTangentObject.computeNumericalTangent
            figName = sprintf('%s_%s_%s_%s_numTang_%s.fig', ...
                setupObject.saveObject.fileName, ...
                continuumObject.elementDisplacementType, ...
                continuumObject.materialObject.name, ...
                setupObject.integrator, ...
                plotTypes{i});
        else
            figName = sprintf('%s_%s_%s_%s_%s.fig', ...
                setupObject.saveObject.fileName, ...
                continuumObject.elementDisplacementType, ...
                continuumObject.materialObject.name, ...
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
        % case 'totalEnergy'
        %     plot(timeVector, internalEnergy, 'color', color_scheme{1});
        %     hold on;
        %     plot(timeVector, kineticEnergy, 'color', color_scheme{5});
        %     plot(timeVector, totalEnergy, 'color', color_scheme{7});
        %     legend('$E_{\mathrm{int}}$', '$E_{\mathrm{kin}}$', '$E_{\mathrm{tot}}$','Interpreter', 'latex');
        %     title(['internal energy, kinetic energy and total energy with $\beta = $', ...
        %         num2str(continuumObject.materialObject.beta), ...
        %         ', $\Delta T$ = ', num2str(timeStepSize), ...
        %         ', $T$ = ', num2str(setupObject.totalTime)], ...
        %    'Interpreter', 'latex');
        case 'totalEnergy'
            plot(timeVector, internalEnergy, 'color', color_scheme{1});
            hold on;
            plot(timeVector, kineticEnergy, 'color', color_scheme{5});
            plot(timeVector, totalEnergy, 'color', color_scheme{7});
            legend('$E_{\mathrm{int}}$', '$E_{\mathrm{kin}}$', '$E_{\mathrm{tot}}$','Interpreter', 'latex');
            title(['internal energy, kinetic energy and total energy', ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], ...
           'Interpreter', 'latex');
        case 'totalEnergyDiff'
            plot(timeVector(tStartDiff + 3:end), totalEnergyDiff, 'color', color_scheme{3});
            legend('$ \Delta E_{\mathrm{tot}}$','Interpreter', 'latex');
            title(['Total energy diff', ... % with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'linearMomentum'
            plot(timeVector, linearMomentum(:,1),'color', color_scheme{1});
            hold on;
            plot(timeVector, linearMomentum(:,2),'color', color_scheme{5});
            plot(timeVector, linearMomentum(:,3),'color', color_scheme{7});
            title(['linear momentum', ...% with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'totalLinearMomentum'
            plot(timeVector, totalLinearMomentum,'color', color_scheme{5});
            title(['total linear momentum', ...% with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'angularMomentum'
            plot(timeVector, angularMomentum(:,1),'color', color_scheme{1});
            hold on;
            plot(timeVector, angularMomentum(:,2),'color', color_scheme{5});
            plot(timeVector, angularMomentum(:,3),'color', color_scheme{7});
            title(['angular momentum', ...% with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'totalAngularMomentum'
            plot(timeVector, totalAngularMomentum, 'color', color_scheme{7});
            title(['total angular momentum', ...% with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'totalAngularMomentumDiff'
            plot(timeVector(tStartDiff + 3:end), totalAngularMomentumDiff, 'color', color_scheme{3});
            legend('$ \Delta J$','Interpreter', 'latex');
            title(['Total angular momentum diff', ... % with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'totalEntropy'
            plot(timeVector, totalEntropy, 'color', color_scheme{1});
            legend('$ S_{\mathrm{tot}}$','Interpreter', 'latex');
            title(['total entropy with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
        case 'totalEntropyDiff'
            plot(timeVector(2:end), totalEntropyDiff(2:end), 'color', color_scheme{3});
            legend('$ \Delta S_{\mathrm{tot}}$','Interpreter', 'latex');
            title(['total entropy diff with $\beta = $', num2str(continuumObject.materialObject.beta), ...
                ', $\Delta T$ = ', num2str(timeStepSize), ...
                ', $T$ = ', num2str(setupObject.totalTime)], 'Interpreter', 'latex');
    end

    % save the plot
    savefig(figPath);
    extraOptions = 'every axis plot/.append style={line width=\figLineWidth}';
    matlab2tikz('height', '\figH', 'width', '\figW','filename', tikzPath, 'showInfo', false, ...
        'extraAxisOptions', {'title style={font=\fontsize{\figFontSize}{\figFontSize}\selectfont}', ...
        'legend style={font=\fontsize{\figFontSize}{\figFontSize}\selectfont}', extraOptions});
    if i < length(plotTypes)
        figure; % create new figure for next plot)
    end
end
