if setupObject.plotFlag
    %% create figure
    if ~exist('figPlot','var')
        scrsz = get(groot,'screensize');
        figPlot = figure('outerposition',[1920,1080,1920,1080]);
    end
    %% plot
    clf
    hold on
    for index1 = 1:dofObject.numberOfContinuumObjects
        if isa(dofObject.listContinuumObjects{index1},'solidClass') || isa(dofObject.listContinuumObjects{index1},'solidThermoClass')
            plot(dofObject.listContinuumObjects{index1})    
        end
    end
    view(3);
    axis([-0.5 4.5 -0.5 4.5 -0.5 4.5]);
    drawnow
end