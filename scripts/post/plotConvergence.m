%% setup
% fileName = 'rotatingXElectroThermoMechanicsDT';
fileName = 'LShapeElectroThermoConvergence';
simulationDuration = 0.1;
variableTimeStepSizeVector = [0.025;0.01;0.005;0.0025;0.001;0.0005]; %[0.1;0.05;0.025;0.01;0.005;0.0025;0.001;0.0005];

%% Compute L2 norm
% load reference solution
lastTimeStep = simulationDuration / min(variableTimeStepSizeVector);
load(strcat(fileName, num2str(min(variableTimeStepSizeVector)), '.mat'), strcat('setupObject', num2str(lastTimeStep)), strcat('dofObject', num2str(lastTimeStep)));
eval(strcat('setupObjectReference = setupObject', num2str(lastTimeStep), ';'));
eval(strcat('dofObjectReference = dofObject', num2str(lastTimeStep), ';'));

continuumObjectReference = dofObjectReference.listContinuumObjects{1};

newTimeStepSizeVector = variableTimeStepSizeVector;
newTimeStepSizeVector(newTimeStepSizeVector == min(variableTimeStepSizeVector)) = [];

errorL2 = struct();
for kk = 1:length(newTimeStepSizeVector)
    lastTimeStep = simulationDuration / newTimeStepSizeVector(kk);
    load(strcat(fileName, num2str(newTimeStepSizeVector(kk)), '.mat'), strcat('setupObject', num2str(lastTimeStep)), strcat('dofObject', num2str(lastTimeStep)));
    eval(strcat('setupObject = setupObject', num2str(lastTimeStep), ';'));
    eval(strcat('dofObject = dofObject', num2str(lastTimeStep), ';'));
    continuumObject = dofObject.listContinuumObjects{1};
    errorL2 = errorCalculationRoutine(errorL2, continuumObject, continuumObjectReference, kk);

    clear dofObject setupObject;
end
variableNames = fieldnames(errorL2);
numberOfVariables = length(variableNames);

%% Plotting
close all
scrsz = get(groot, 'screensize');
bplot = scrsz(3) * 0.2;
hplot = scrsz(4) * 0.4;
xplot = scrsz(3) * [0.01, 0.32, 0.63];
yplot = scrsz(4) * [0.5, 0.1];
figPlot = figure('outerposition', [xplot(3), yplot(2), bplot, hplot]);
hold on
% symbols = ['*'; '+'; 'o'; 'x'];
AT = newTimeStepSizeVector;
slopes = zeros(numberOfVariables, 1);
for markerCounter = 1:numberOfVariables
    plot(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))), '-x')
    % plot(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))), strcat('-', symbols(markerCounter)));
    %     loglog(AT,vertcat(errorL2.(variableNames{markerCounter})),strcat('-',symbols(markerCounter)));
    yp = polyfit(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))), 1);
    slopes(markerCounter) = yp(1);
end
axis tight;
grid on;
% set(gca,'xscale','log','yscale','log')
title('Log_{10} of L2-Norm of error');
xlabel('\Delta t');
ylabel('error');
% legendEntries = cell(numberOfVariables, 1);
legendEntries = variableNames;
for ii=1:numberOfVariables
    switch variableNames{ii}
        case 'electroPhi'
            legendEntries{ii} = '\Phi';
        case 'theta'
            legendEntries{ii} = '\theta';
        case 'lambdaC'
            legendEntries{ii} = '\Lambda^C';
        case 'lambdaG'
            legendEntries{ii} = '\Lambda^G';
        case 'lambdac'
            legendEntries{ii} = '\Lambda^c';
    end
end
% legend('x', '\Phi', '\theta', 'D') %,'location','west')
legend(legendEntries);
set(gca, 'fontsize', 16)
% a1 = 1;
% ae = size(errorL2, 1);
% propGrad = (log10(errorL2(a1, :)) - log10(errorL2(ae, :))) ./ (log10(1./AT(a1)) - log10(1./AT(ae)));

disp('Slopes:');
for ii=1:numberOfVariables
    disp(strcat(variableNames{ii}, ':  ',num2str(slopes(ii))));
end
