if setupObject.plotFlag
    %% create figure
    if ~exist('figPlot','var')
        scrsz = get(groot,'screensize');
        figPlot = figure('outerposition',[1920,1080,1920,1080]);
    end
    %% plot
    clf
    if ismember('solidObjects',evalin('base','who'))
        plot(evalin('base','solidObjects'))    
    elseif ismember('solidThermoObjects',evalin('base','who'))
        plot(evalin('base','solidThermoObjects'))
    end
    view(3);
    axis([-0.5 4.5 -0.5 4.5 -0.5 4.5]);
    drawnow
end