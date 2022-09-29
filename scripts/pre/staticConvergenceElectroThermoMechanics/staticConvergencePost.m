close all
% get all discretisations
flagPlot = true;

errorL2 = zeros(size(All_Elem,2),4);
%% Loop over all elements
for ii = 1:numel(choiceElement)
    %% Loop over all discretisations
    for j = 1:size(All_Elem,2)
        %         load(strcat('staticConvergence', choiceElement{ii}, '.mat'))
        %         load(strcat('staticConvergence', choiceElement{ii}, '.mat'),'journalSolidElectroThermoObject','analyticalVariablesCell')
        out = error_calculation_routine(solidElectroThermoObjectCell{ii,j},analyticalVariablesCell{ii,j});
        deltaL2 = sqrt(sum(vertcat(out.delta)));
        anaL2 = sqrt(sum(vertcat(out.ana)));
        %Displacement L2-Error
        errorL2(j,1) = sqrt(sum(vertcat(out.delta)))/sqrt(sum(vertcat(out.ana)));
        errorL2(j,2) = sqrt(sum(vertcat(out.delta_phi)))/sqrt(sum(vertcat(out.phi_analytical)));
        errorL2(j,3) = sqrt(sum(vertcat(out.delta_theta)))/sqrt(sum(vertcat(out.theta_analytical)));
        errorL2(j,4) = sqrt(sum(vertcat(out.delta_D)))/sqrt(sum(vertcat(out.D_analytical)));
        errorL2(j,5) = sqrt(sum(vertcat(out.delta_C)))/sqrt(sum(vertcat(out.C_analytical)));
        errorL2(j,6) = sqrt(sum(vertcat(out.delta_G)))/sqrt(sum(vertcat(out.G_analytical)));
        errorL2(j,7) = sqrt(sum(vertcat(out.delta_c)))/sqrt(sum(vertcat(out.c_analytical)));
        errorL2(j,8) = sqrt(sum(vertcat(out.delta_lambdaC)))/sqrt(sum(vertcat(out.lambdaC_analytical)));
        errorL2(j,9) = sqrt(sum(vertcat(out.delta_lambdaG)))/sqrt(sum(vertcat(out.lambdaG_analytical)));
        errorL2(j,10) = sqrt(sum(vertcat(out.delta_lambdac)))/sqrt(sum(vertcat(out.lambdac_analytical)));
    end
    %% Plotting
    scrsz = get(groot,'screensize');
    bplot = scrsz(3)*0.2;
    hplot = scrsz(4)*0.4;
    xplot = scrsz(3)*[0.01,0.32,0.63];
    yplot = scrsz(4)*[0.5,0.1];
    figPlot = figure('outerposition',[xplot(3),yplot(2),bplot,hplot]);
    hold on
    symbols = ['*';'+';'d';' ';'s';'o';'x';'.';'-';'v'];
    for markerCounter = 1:numel(symbols)
        plot(log10(1./All_Elem),log10(errorL2(:,markerCounter)),strcat('-',symbols(markerCounter)))
    end
    axis tight;
    axis([-0.925 -0.6 -8 0])
    grid on;
    %set(gca,'xscale','log','yscale','log')
    title('Log_{10} of L2-Norm of error');
    xlabel('Log_{10}(h)');
    ylabel('Log_{10} of L2-Norm of error');
    legend('x','\Phi','\theta','D','C','G','c','\Lambda_C','\Lambda_G','\Lambda_c')%,'location','west')
    set(gca,'fontsize',16)
    a1 = 1;
    ae = size(errorL2,1);
    propGrad = (log10(errorL2(a1,:))-log10(errorL2(ae,:)))./(log10(1./All_Elem(a1))-log10(1./All_Elem(ae)))
    try
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
    catch
    end

    % Exporting
    filename = fullfile(strcat('hconvergence',choiceElement{ii}));
    if flagPlot
        print(figPlot,filename,'-dpng')
        print(figPlot,filename,'-depsc')
        savefig(figPlot,filename)
        matlab2tikz([filename,'.tikz'],'width','\figW','height','\figH')
    end
    clf(figPlot)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     continue
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %% Plot initial and final configuration
    solidElectroThermoObject = solidElectroThermoObjectCell{ii,j};
    if 1
        setupObject.plotObject.time = 'R';
        %
        filename = fullfile(strcat('initialConfiguration',choiceElement{ii}));
        setupObject.plotObject.postPlotType = 'zero';
        plot(solidElectroThermoObject,setupObject)
        axis([-0.5 1.5 -0.5 1.5 -0.5 1.5])
        axis equal
        axis tight
        axis off
        view(3)
        colorbar off
        if flagPlot
            print(figPlot,filename,'-depsc')
            print(figPlot,filename,'-dpng')
            savefig(figPlot,filename)
        end
        clf(figPlot)
    end

    setupObject.plotObject.time = 'N1';
    plotVarVector{1} = 'phi';
    plotVarVector{2} = 'temp';
    plotVarVector{3} = 'stress';
    for jjj = 1:2
        if jjj == 1
            flagPlotColorbar = false;
        elseif jjj == 2
            flagPlotColorbar = true;
        end
        for jj = 1:1:numel(plotVarVector)
            setupObject.plotObject.postPlotType = plotVarVector{jj};
            plot(solidElectroThermoObject,setupObject)
            axis([-0.5 1.5 -0.5 1.5 -0.5 1.5])
            axis equal
            axis tight
            view(3)
            colorbar off
            if flagPlotColorbar
                clf(figPlot)
                filename = fullfile(strcat(plotVarVector{jj},choiceElement{ii},'Colorbar'));
            end
            c = colorbar('Location','south');
            axis off
            if strcmpi(plotVarVector{jj},'stress')
                colormap(parula)
                set(gca,'clim',[0 2e4]);
            elseif strcmpi(plotVarVector{jj},'phi')
                colormap(jet)
                set(gca,'clim',[0 100]);
            elseif strcmpi(plotVarVector{jj},'temp')
                colormap(flipud(autumn))
                set(gca,'clim',[293 303]);
            end
            if ~flagPlotColorbar
                colorbar off
                filename = fullfile(strcat(plotVarVector{jj},choiceElement{ii}));
            end
            set(gca,'fontsize',16)
            if flagPlot
                matlab2tikz([filename,'.tikz'],'width','\figW','height','\figH')
                print(figPlot,filename,'-depsc')
                print(figPlot,filename,'-dpng')
                savefig(figPlot,filename)
            end
            clf(figPlot)
        end
    end
end