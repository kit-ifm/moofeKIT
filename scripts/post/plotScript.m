if setupObject.plotObject.flag
    if setupObject.plotObject.everyNewtonStep
        %% plot
        if ~exist('figurePlot','var')
            if ~setupObject.plotObject.docked
                scrsz = get(groot,'screensize');
                figurePlot = figure('outerposition',[scrsz(3)*.4,scrsz(4)*.1,scrsz(3)*.6,scrsz(4)*.85]);
            else
                figurePlot = figure();
            end
        else
            cla(figureAxes)
        end
        plotObject = setupObject.plotObject;
        %formatting of axis
        if ~exist('timeStep','var') || timeStep==1
            plotObject.findXminXmax(dofObject.listContinuumObjects{1})  % FIXME listContinuumObjects{1}: other objects?
            %scale axes
            xlim([plotObject.xmin(1)-plotObject.deltaXY*plotObject.border, plotObject.xmax(1)+plotObject.deltaXY*plotObject.border]);
            ylim([plotObject.xmin(2)-plotObject.deltaXY*plotObject.border, plotObject.xmax(2)+plotObject.deltaXY*plotObject.border]);
            if dofObject.listContinuumObjects{1}.dimension == 3
                zlim([plotObject.xmin(3)-plotObject.deltaXY*plotObject.border, plotObject.xmax(3)+plotObject.deltaXY*plotObject.border]);
            end
            axis equal
            figureAxes = gca;
            % ensure that axis limits stay the same
            figureAxes.XLimMode = 'manual';
            figureAxes.YLimMode = 'manual';
            if dofObject.listContinuumObjects{1}.dimension == 3
                figureAxes.ZLimMode = 'manual';
            end
            hold on
            if ~plotObject.axis
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                if dofObject.listContinuumObjects{1}.dimension == 3
                    set(gca,'zticklabel',[])
                    set(gca,'ztick',[])
                end
            end
            view(plotObject.view)
        end
        %% plot objects
        for index1 = 1:dofObject.numberOfContinuumObjects
            % solids
            if isa(dofObject.listContinuumObjects{index1}, 'solidSuperClass')
                plot(dofObject.listContinuumObjects{index1},setupObject);
            end
            % boundarys 2D
            if dofObject.listContinuumObjects{index1}.dimension == 2
                if isa(dofObject.listContinuumObjects{index1},'dirichletClass') || isa(dofObject.listContinuumObjects{index1},'neumannClass')
                    plot(dofObject.listContinuumObjects{index1},setupObject);
                end
            end
        end
        if ~isempty(setupObject.plotObject.colorBarLimits)
            caxis(setupObject.plotObject.colorBarLimits)
        end
        drawnow
    end
end
% TODO: configuration to export multiple vtk files
if setupObject.plotObject.makeMovie
    generatedFilesFolderInfo = what('generatedFiles');
    filename = strcat(generatedFilesFolderInfo.path, filesep, 'movieFiles', filesep, setupObject.saveObject.fileName, 'Movie', num2str(timeStep), '.vtk');
    for index1 = 1:dofObject.numberOfContinuumObjects
        obj = dofObject.listContinuumObjects{index1};
        if isa(obj, 'solidSuperClass')
            if isa(obj,'solidElectroClass')
                VTKPlot(filename, 'unstructured_grid', obj.qN1(:,1:3), obj.edof, 'scalars', 'electricPotential', obj.qN1(:,4));
            elseif isa(obj,'solidThermoClass')
                VTKPlot(filename, 'unstructured_grid', obj.qN1(:,1:3), obj.edof, 'scalars', 'temperature', obj.qN1(:,4));
            elseif isa(obj,'solidElectroThermoClass')
%                 VTKPlot(filename,'unstructured_grid',obj.qN1(:,1:3),obj.meshObject.edof,'scalars','temperature',obj.qN1(:,5))
                VTKPlot(filename, 'unstructured_grid', obj.qN1(:,1:3), obj.meshObject.edof, 'scalars', 'electricPotential', obj.qN1(:,4));
            else
                error(['Movie maker is not implemented for objects of class', ' ', class(obj), '!']);
            end
        end
    end
end