%% Script for plotting the total energy
clear;
clc;

% general settings
fileNames = {'rotatingXElectroThermoMechanicsMidpoint'};%, 'rotatingXMidpoint'};
plotLegend = true;
legendNames = {'MP'};
legendLocation = 'east';
plotVerticalLines = false;
plotVerticalLinesAtTime = [0.4, 0.5, 0.9];

% settings total energy plot
exportAsTikz = false;

% settings total energy difference plot
plotTotalEnergyDifferences = true;
plotTotalEnergyDifferencesForFile = 1;
plotOnlyIncrementalEnergyChanges = 0;
exportTotalEnergyDifferencesAsTikz = true;

%% plot total energy

% load data
numberOfPlots = length(fileNames);
setupObjects = cell(numberOfPlots, 1);
dofObjects = cell(numberOfPlots, 1);
for ii=1:numberOfPlots
    [setupObjectCell, dofObjectCell] = loadObjectsFromMatFile(fileNames{ii}, 'last');
    setupObjects{ii} = setupObjectCell{1};
    dofObjects{ii} = dofObjectCell{1};
end

% plot data
timeVectors = cell(numberOfPlots, 1);
totalEnergies = cell(numberOfPlots, 1);
for ii=1:numberOfPlots
    setupObject = setupObjects{ii};
    dofObject = dofObjects{ii};
    timeVector = getTime(dofObject.postDataObject, setupObject);
    kineticEnergy = getKineticEnergy(dofObject.postDataObject, setupObject);
    internalEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'internalEnergy');
    externalEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'externalEnergy');
    totalEnergy = kineticEnergy + internalEnergy;
    if ii == 1
        figure;
    end
    plot(timeVector, totalEnergy);
    timeVectors{ii} = timeVector;
    totalEnergies{ii} = totalEnergy;
    if ii == 1
        xlim([timeVector(1), timeVector(end)]);
        xlabel('$t$ [s]');
        ylabel('$E_{n+1}$ [J]');
        hold on;
    end
end

% plot vertical lines
if plotVerticalLines
    yLimit = ylim;
    for ii = 1:length(plotVerticalLinesAtTime)
        plot(plotVerticalLinesAtTime(ii)*ones(2,1), [yLimit(1); yLimit(2)], '--k');
    end
end

% plot legend
if plotLegend
    legend(legendNames, 'Location', legendLocation);
end

%% export total energy as tikz
if exportAsTikz
    matlab2tikz(strcat(filename, 'TotalEnergy', '.tikz'), 'width', '\fwidth', 'height', '\fheight');
end

%% plot total energy differences
if plotTotalEnergyDifferences
    setupObject = setupObjects{plotTotalEnergyDifferencesForFile};
    dofObject = dofObjects{plotTotalEnergyDifferencesForFile};
    timeVector = timeVectors{plotTotalEnergyDifferencesForFile};
    totalEnergy = totalEnergies{plotTotalEnergyDifferencesForFile};

    totalEnergyBefore = totalEnergy(1:end-1);
    totalEnergyDifference = totalEnergy(2:end) - totalEnergyBefore;
    timeVector = timeVector(1:end-1);

    if plotOnlyIncrementalEnergyChanges
        newtonTolerance = setupObject.newton.tolerance;
        for ii = 1:size(totalEnergyDifference, 1)
            if abs(totalEnergyDifference(ii)) > newtonTolerance
                totalEnergyDifference(ii) = NaN;
            end
        end
    end
    figure;
    plot(timeVector, totalEnergyDifference);
    yLimit = 2*max(totalEnergyDifference);

    if isnan(yLimit)
        error("Can't plot total energy differences. The difference of total energy between one time step and another is too large!")
    end

    xlim([timeVector(1), timeVector(end)]);
    ylim([-yLimit, yLimit]);
    xlabel('$t$ [s]');
    ylabel('$E_{n+1}-E_n$ [J]');

    % plot horizontal line at y=0
    hold on;
    plot(timeVector, zeros(size(timeVector, 1), 1), '--k');

    % plot vertical lines
    if plotVerticalLines
        for ii = 1:length(plotVerticalLinesAtTime)
            plot(plotVerticalLinesAtTime(ii)*ones(2,1), [-yLimit; yLimit], '--k');
        end
    end
end