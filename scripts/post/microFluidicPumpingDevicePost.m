%% postprocessing
% load data
loadingVector = [1, 60, 220, 1000];
[setupObjectCell, dofObjectCell] = loadObjectsFromMatFile('microFluidicPumpingDevicePre',loadingVector);

setupObject = setupObjectCell{1};
dofObject = dofObjectCell{1};
solidElectroThermoObject = dofObject.listContinuumObjects{1};
if 0
    %% plots
    % phi
    figure;
    dt = setupObject.timeStepSize;
    totalTime = setupObject.totalTime;
    plot([0:dt:totalTime],dirichletObject1a.timeFunction([0:dt:totalTime]),'linewidth',2);
    hold on;
    plot([0:dt:totalTime],dirichletObject1b.timeFunction([0:dt:totalTime]),'linewidth',2);
    legend('phi2','phi3')
    xlabel('t [s]');
    ylabel('phi [V]');
    matlab2tikz(['phi','.tikz'],'width','\fwidth','height','\fheight')
    % Q
    figure;
    plot([0:dt:totalTime],neumannObject1.energy*neumannObject1.timeFunction([0:dt:totalTime]),'k','linewidth',2);
    xlabel('t [s]');
    ylabel('Q ');
    matlab2tikz(['Q','.tikz'],'width','\fwidth','height','\fheight')
    % assign data
    timeVector = getTime(dofObject.postDataObject,setupObject);
    kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
    internalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'internalEnergy');
    DE = getEnergy(dofObject.postDataObject,dofObject,setupObject,'DE');
    TS = getEnergy(dofObject.postDataObject,dofObject,setupObject,'TS');
    externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
    [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
    [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
    % energy
    figure;
    totalEnergy = kineticEnergy + internalEnergy;
    plot(timeVector,totalEnergy);
    xlabel('$t$ [s]');
    ylabel('$E_{n+1}$ [J]');
    matlab2tikz(['energy','.tikz'],'width','\fwidth','height','\fheight')
    % energy difference
    figure; hold on;
    timeVectorDiff = timeVector(2:end);
    newtonTolerance = setupObject.newton.tolerance;
    totalEnergyN = totalEnergy(1:end-1);
    totalEnergyDiff = totalEnergy(2:end) - totalEnergyN;
    for ii = 1:(numel(bcTimeEnd)-2)/2
        diffIndex = int32(bcTimeEnd(2*ii)/dt+1:bcTimeEnd(2*ii+1)/dt);
        plot(timeVectorDiff(diffIndex), totalEnergyDiff(diffIndex), 'k','linewidth',2);
    end
    plotHelp(1) = plot([dt,timeVectorDiff(end)],[0,0],'k--');
    for jj = 1:numel(bcTimeEnd)-2
        plotHelp(jj+1) = plot([bcTimeEnd(jj+1),bcTimeEnd(jj+1)],newtonTolerance*[-1,1],'k-.');
    end
    box on
    xlabel('$t$ [s]');
    ylabel('$E_{n+1}-E_{n}$ [J]');
    set(gca,'fontsize',16)
    axis([0 dt*setupObject.totalTimeSteps -newtonTolerance newtonTolerance])
    matlab2tikz(['energyDifference','.tikz'],'width','\fwidth','height','\fheight')
    % linear momentum p = m*v [N*s]
    % figure;
    % plot(timeVector,totalLinearMomentum);
    % angular momentum j = rxp [J*s]
    % figure;
    % plot(timeVector,totalAngularMomentum);
    % xlabel('$t$ [s]');
    % ylabel('$J_{n+1}$ [Js]');
    % % angular momentum difference
    % figure;
    % totalAngularMomentumN = totalAngularMomentum(1:end-1);
    % totalAngularMomentumDiff = totalAngularMomentum(2:end) - totalAngularMomentumN;
    % plot(timeVector(2:end), totalAngularMomentumDiff);
    % xlabel('$t$ [s]');
    % ylabel('$J_{n+1}-J_{n}$ [J]');
end

%% snapshots
% initial postion
scrsz = get(groot,'screensize');
bplot = scrsz(3)*0.2;
hplot = scrsz(4)*0.4;
xplot = scrsz(3)*[0.01,0.32,0.63];
yplot = scrsz(4)*[0.5,0.1];
figPlot = figure('outerposition',[xplot(3),yplot(2),bplot,hplot]);
setupObject.plotObject.time = 'R';
plot(solidElectroThermoObject,setupObject)
axis equal
axis tight
axis off
view(-8,13)
axis([0 0.55 0 0.5 0 0.2])
colorbar off
print(figPlot,'0','-depsc')
% print(figPlot,'0','-dpng')
% savefig(figPlot,'0')
clf(figPlot)

colorBarCell = {'temp', 'stress', 'phi', 'D1';...
                 [293 346], [0 40], [-phi0 phi0], [-4 4]*1e-10;...
                cellstr('flipud(autumn)'), 'parula', 'jet', 'spring'};
flagInitial = false;

for ii = 1:size(setupObjectCell,1)
    for jj = 1:size(colorBarCell,2)
        for jjj = 1
            switch jjj 
                case 1
                    flagPlotColorbar = false;
                case 2
                    flagPlotColorbar = true;
            end
            setupObject = setupObjectCell{ii};
            dofObject = dofObjectCell{ii};
            solidElectroThermoObject = dofObject.listContinuumObjects{1};
            setupObject.plotObject.postPlotType = colorBarCell{1,jj};
            if ~flagInitial
                setupObject.plotObject.time = 'N1';
            else
                setupObject.plotObject.time = 'R';
            end
            if strcmpi(setupObject.plotObject.time,'R')
                saveName = strcat(colorBarCell{1,jj},setupObject.plotObject.time);
            else
                saveName = strcat(colorBarCell{1,jj},num2str(loadingVector(ii)));
            end
            %% plot half pump
            plot(solidElectroThermoObject,setupObject)
            axis equal
            axis tight
            axis off
            axis([-0.55 0.55 -0.55 0.55 -0.25 0.25])
            % view(3)
            view(-8,13)
            hold on
            %
            if ~flagInitial
                for kk = 1:3
                    %                 solidElectroThermoObjectArray(kk) = copyContinuumObject(dofObject,solidElectroThermoObject);
                    solidElectroThermoObjectArray(kk) = copy(solidElectroThermoObject);
                    switch kk
                        case 1
                            solidElectroThermoObjectArray(kk).reflection(3)
                        case 2
                            solidElectroThermoObjectArray(kk).rotate(2,pi)
                        case 3
                            solidElectroThermoObjectArray(kk).reflection(1)
                    end
                    plot(solidElectroThermoObjectArray(kk),setupObject)
                end
            else
                for kk = 1:7
                    %                 solidElectroThermoObjectArray(kk) = copyContinuumObject(dofObject,solidElectroThermoObject);
                    solidElectroThermoObjectArray(kk) = copy(solidElectroThermoObject);
                    switch kk
                        case 1
                            solidElectroThermoObjectArray(kk).reflection(3)
                        case 2
                            solidElectroThermoObjectArray(kk).rotate(2,pi)
                        case 3
                            solidElectroThermoObjectArray(kk).reflection(1)
                        case 4
                            solidElectroThermoObjectArray(kk).reflection(2)
                        case 5
                            solidElectroThermoObjectArray(kk).rotate(1,pi)
                        case 6
                            solidElectroThermoObjectArray(kk).rotate(3,pi)
                        case 7
                            solidElectroThermoObjectArray(kk).rotate(1,pi)
                            solidElectroThermoObjectArray(kk).reflection(1)
                    end
                    plot(solidElectroThermoObjectArray(kk),setupObject)
                end
                saveName = strcat(saveName,'Initial');
            end
            if flagPlotColorbar
                clf(figPlot)
                saveName = strcat(saveName,'Colorbar');
            end            
            c = colorbar('Location','south');
            axis off
            if jj==1
                colormap(flipud(autumn))
            else
                colormap(colorBarCell{3,jj})
            end
            caxis(colorBarCell{2,jj});
            set(gca,'fontsize',16)
            if ~flagPlotColorbar
                colorbar off
            end
%             print(figPlot,saveName,'-depsc');
            print(figPlot,saveName,'-dpng');
            if flagPlotColorbar
                matlab2tikz([saveName,'.tikz'],'width','\figW','height','\figH')
            end
            clf(figPlot);
        end
    end
end