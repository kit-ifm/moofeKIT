close all
clearvars figPlot
% get all discretisations 
errorL2 = zeros(size(All_Elem,2),4);
%% Loop over all discretisations
for j = 1:size(All_Elem,2)
    out = error_calculation_routine(journalSolidElectroThermoObject{j},analyticalVariablesCell{j});
    deltaL2 = sqrt(sum(vertcat(out.delta)));
    anaL2 = sqrt(sum(vertcat(out.ana)));
    %Displacement L2-Error
    errorL2(j,1) = sqrt(sum(vertcat(out.delta)))/sqrt(sum(vertcat(out.ana)));
    errorL2(j,2) = sqrt(sum(vertcat(out.delta_C)))/sqrt(sum(vertcat(out.C_analytical)));
    errorL2(j,3) = sqrt(sum(vertcat(out.delta_G)))/sqrt(sum(vertcat(out.G_analytical)));
    errorL2(j,4) = sqrt(sum(vertcat(out.delta_c)))/sqrt(sum(vertcat(out.c_analytical)));
    errorL2(j,5) = sqrt(sum(vertcat(out.delta_lambdaC)))/sqrt(sum(vertcat(out.lambdaC_analytical)));
    errorL2(j,6) = sqrt(sum(vertcat(out.delta_lambdaG)))/sqrt(sum(vertcat(out.lambdaG_analytical)));
    errorL2(j,7) = sqrt(sum(vertcat(out.delta_lambdac)))/sqrt(sum(vertcat(out.lambdac_analytical)));
end

%% Plotting
close all
scrsz = get(groot,'screensize');
bplot = scrsz(3)*0.2;
hplot = scrsz(4)*0.4;
xplot = scrsz(3)*[0.01,0.32,0.63];
yplot = scrsz(4)*[0.5,0.1];
figPlot = figure('outerposition',[xplot(3),yplot(2),bplot,hplot]);
hold on
symbols = ['*';'+';'d';' ';'s';'o';'x'];
for markerCounter = 1:numel(symbols)
    plot(log10(1./All_Elem),log10(errorL2(:,markerCounter)),strcat('-',symbols(markerCounter)))
end
axis tight;
grid on;
%set(gca,'xscale','log','yscale','log')
title('Log_{10} of L2-Norm of error');
xlabel('Log_{10}(h)');
ylabel('Log_{10} of L2-Norm of error');
legend('x','C','G','c','\Lambda_C','\Lambda_G','\Lambda_c')%,'location','west')
set(gca,'fontsize',16)
a1 = 1;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))



a1 = 2;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))




a1 = 3;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))

a1 = 4;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))

a1 = 1;
ae = size(errorL2,1)-1;
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))



a1 = 6;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))
a1 = 7;
ae = size(errorL2,1);
propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))
% Exporting
filename = fullfile('hconvergence');
print(figPlot,filename,'-dpng')
print(figPlot,filename,'-depsc')
savefig(figPlot,filename)
matlab2tikz([filename,'.tikz'],'width','\figW','height','\figH')
clf(figPlot)

%% Plot initial and final configuration
plotTime = 'REF';
flagPlotColorbar = false;
% 
filename = fullfile(savepath,'initialConfiguration');
% plot(solidElectroThermoSystem,'mesh','time','REF','view',solidElectroThermoSystem.DIM)
staticConvergencePlotTool
print(figPlot,filename,'-depsc')
print(figPlot,filename,'-dpng')
savefig(figPlot,filename)
clf(figPlot)

plotTime = 'N1';
flagPlotColorbar = false;
plotVarVector{1} = 'phi';
plotVarVector{2} = 'temp';
plotVarVector{3} = 'stress';
for ii = 1:1:numel(plotVarVector)
    plotVar = plotVarVector{ii};
    staticConvergencePlotTool
    if flagPlotColorbar
        filename = fullfile(savepath,strcat(plotVarVector{ii},'Colorbar'));
    else
        filename = fullfile(savepath,plotVarVector{ii});
    end
    if flagPlotColorbar
        matlab2tikz([filename,'.tikz'],'width','\figW','height','\figH')
    else
        colorbar off
    end
    print(figPlot,filename,'-depsc')
    print(figPlot,filename,'-dpng')
    savefig(figPlot,filename)
    clf(figPlot)
end
