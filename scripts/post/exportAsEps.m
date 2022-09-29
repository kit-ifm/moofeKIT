%% Script to create .eps files from .mat files
clear;
clc;

filename = 'rotatingXElectroThermoMechanics';
numberOfTimeSteps = 30;
createEpsForTimeSteps = 1:30;
plotMode = 'onlyContinuumObject';

figureView = [-41.33, 24.725];
figurePosition = [660,615,580,363];
axisValues = [-0.03 0.03 -0.03 0.03 0 0.015];

%% determine values for plot
minX = 1e12;
maxX = -1e12;
minY = 1e12;
maxY = -1e12;
minZ = 1e12;
maxZ = -1e12;
colorBarMin = 1e12;
colorBarMax = -1e12;
listObjects = cell(size(createEpsForTimeSteps));
for ii = 1:numberOfTimeSteps
    if ismember(ii, createEpsForTimeSteps)
        load(strcat(filename, num2str(ii), '.mat'));
        continuumObject = dofObject.listContinuumObjects{1};
        listObjects{ii} = {setupObject, dofObject};
        
        dimension = continuumObject.dimension;
        
        % axis values
        temporaryMinX = min(continuumObject.qN1(:,1));
        temporaryMaxX = max(continuumObject.qN1(:,1));
        temporaryMinY = min(continuumObject.qN1(:,2));
        temporaryMaxY = max(continuumObject.qN1(:,2));
        if dimension == 3
            temporaryMinZ = min(continuumObject.qN1(:,3));
            temporaryMaxZ = max(continuumObject.qN1(:,3));
        end
        
        if temporaryMinX < minX; minX = temporaryMinX; end
        if temporaryMaxX > maxX; maxX = temporaryMaxX; end
        if temporaryMinY < minY; minY = temporaryMinY; end
        if temporaryMaxY > maxY; maxY = temporaryMaxY; end
        if dimension == 3
            if temporaryMinZ < minZ; minZ = temporaryMinZ; end
            if temporaryMaxZ > maxZ; maxZ = temporaryMaxZ; end
        end
        
        % color bar values
        setupObject.plotObject.postPlotType = 'phi';
        if strcmpi(setupObject.plotObject.postPlotType, 'stress')
            % TODO
        elseif strcmpi(setupObject.plotObject.postPlotType, 'temp')
            if isa(continuumObject, 'solidThermoClass')
                temporaryColorBarMin = min(continuumObject.qN1(:,dimension+1));
                temporaryColorBarMax = max(continuumObject.qN1(:,dimension+1));
            elseif isa(continuumObject, 'solidElectroThermoClass')
                temporaryColorBarMin = min(continuumObject.qN1(:,dimension+2));
                temporaryColorBarMax = max(continuumObject.qN1(:,dimension+2));
            end
        elseif strcmpi(setupObject.plotObject.postPlotType, 'phi')
            temporaryColorBarMin = min(continuumObject.qN1(:,dimension+1));
            temporaryColorBarMax = max(continuumObject.qN1(:,dimension+1));
        end
        
        if temporaryColorBarMin < colorBarMin; colorBarMin = temporaryColorBarMin; end
        if temporaryColorBarMax > colorBarMax; colorBarMax = temporaryColorBarMax; end
        
        % TODO: code duplication here...
        % take into account initial configuration
        if ii == 1
            % axis values
            temporaryMinX = min(continuumObject.qR(:,1));
            temporaryMaxX = max(continuumObject.qR(:,1));
            temporaryMinY = min(continuumObject.qR(:,2));
            temporaryMaxY = max(continuumObject.qR(:,2));
            if dimension == 3
                temporaryMinZ = min(continuumObject.qR(:,3));
                temporaryMaxZ = max(continuumObject.qR(:,3));
            end

            if temporaryMinX < minX; minX = temporaryMinX; end
            if temporaryMaxX > maxX; maxX = temporaryMaxX; end
            if temporaryMinY < minY; minY = temporaryMinY; end
            if temporaryMaxY > maxY; maxY = temporaryMaxY; end
            if dimension == 3
                if temporaryMinZ < minZ; minZ = temporaryMinZ; end
                if temporaryMaxZ > maxZ; maxZ = temporaryMaxZ; end
            end

            % color bar values
            if strcmpi(setupObject.plotObject.postPlotType, 'stress')
                % TODO
            elseif strcmpi(setupObject.plotObject.postPlotType, 'temp')
                if isa(continuumObject, 'solidThermoClass')
                    temporaryColorBarMin = min(continuumObject.qR(:,dimension+1));
                    temporaryColorBarMax = max(continuumObject.qR(:,dimension+1));
                elseif isa(continuumObject, 'solidElectroThermoClass')
                    temporaryColorBarMin = min(continuumObject.qR(:,dimension+2));
                    temporaryColorBarMax = max(continuumObject.qR(:,dimension+2));
                end
            elseif strcmpi(setupObject.plotObject.postPlotType, 'phi')
                temporaryColorBarMin = min(continuumObject.qR(:,dimension+1));
                temporaryColorBarMax = max(continuumObject.qR(:,dimension+1));
            end

            if temporaryColorBarMin < colorBarMin; colorBarMin = temporaryColorBarMin; end
            if temporaryColorBarMax > colorBarMax; colorBarMax = temporaryColorBarMax; end
        end
    end
end

%% plot and export
for ii = 1:numberOfTimeSteps
    if ismember(ii, createEpsForTimeSteps)
        % load continuumObject for current time step
        setupObject = listObjects{ii}{1};
        dofObject = listObjects{ii}{2};
        continuumObject = dofObject.listContinuumObjects{1};
        
        setupObject.plotObject.flag = true;
%         setupObject.plotObject.time = 'R';
        setupObject.plotObject.view = figureView;
        plotScript;
        figurePlot = gcf;
        figureAxes = gca;

        % modify the appearance
        if strcmpi(plotMode, 'onlyContinuumObject')
            figurePlot.Color = [1, 1, 1];
            figurePlot.Position = figurePosition;
            colorBar = colorbar;
            colormap(jet)
            caxis([floor(colorBarMin), ceil(colorBarMax)]);
            colorBar.Visible = false;
            if isempty(axisValues)
                axis([minX maxX minY maxY minZ maxZ]);
            else
                axis(axisValues);
            end
            figureAxes.Visible = false;
            figureAxes.Position = [0,0,1,1];
            figureAxes.Clipping = false;
            
            
        end

        % export as .eps
        exportgraphics(figureAxes, strcat(filename, num2str(ii), '.eps'));
        
        % export colorbar as .eps
%         if ii == createEpsForTimeSteps(end)
%             figureAxes.Children.Visible = false;
%             colorBar = colorbar;
%             colorBar.Visible = true;
%             colorBar.Location = 'south';
%             colorBarPosition = colorBar.Position;
%             colorBar.Position = [0.15, 0.5, 0.7, colorBarPosition(4)];
%             colorBarTickLabels = colorBar.TickLabels;
%             for jj = 1:size(colorBarTickLabels, 1)
%                 if jj == 1
%                     colorBarTickLabels{jj} = 'A';
%                 end
%                 if jj ~= 1 && jj ~= (size(colorBarTickLabels, 1)+1)/2 && jj ~= size(colorBarTickLabels, 1)
%                     colorBarTickLabels{jj} = '';
%                 end
%             end
%             colorBar.TickLabels = colorBarTickLabels;
%             
%             exportgraphics(figureAxes, strcat(filename, 'Colorbar.svg'), 'ContentType','vector');
%         end
    end
end
