%% create figure
if ~exist('figPlot','var')
    scrsz = get(groot,'screensize');
    figPlot = figure('outerposition',[1920,1080,1920,1080]);
end
%% plot
clf
plot(solidObject)
view(3);
axis([-0.5 4.5 -0.5 4.5 -0.5 4.5]);
drawnow