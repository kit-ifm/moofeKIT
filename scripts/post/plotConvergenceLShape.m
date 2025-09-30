%% setup
% fileName = 'generatedFiles/matFiles/lShapedisplacementSCpHCGJDiscreteGradientLShapeH1';
fileName = 'generatedFiles/matFiles/lShapemixedSCpHCGJLambdaDiscreteGradientLShapeH1';
simulationDuration = 1;
% variableTimeStepSizeVector = [0.5; 0.25; 0.1; 0.05; 0.025; 0.01; 0.005; 0.0025; 0.001; 0.0005];%; 0.00025; 0.0001; 0.00001];
% variableTimeStepSizeVector = [0.025; 0.01; 0.005; 0.0025; 0.001; 0.0005];%; 0.00025; 0.0001; 0.00001];
variableTimeStepSizeVector = [0.025; 0.01; 0.005; 0.0025; 0.0005];%; 0.00025; 0.0001; 0.00001];
%% Compute L2 norm
% load reference solution
lastTimeStep = simulationDuration / min(variableTimeStepSizeVector);
load(strcat(fileName, autoNum2Str(min(variableTimeStepSizeVector)), '.mat'), strcat('setupObject', num2str(lastTimeStep)), strcat('dofObject', num2str(lastTimeStep)));
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
    errorL2 = errorCalculationRoutineMF(errorL2, continuumObject, continuumObjectReference, kk,setupObject);
    clear dofObject setupObject;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% remove fields %%%%%%%%%%%%
errorL2 = rmfield(errorL2,'lambdaC');
errorL2 = rmfield(errorL2,'lambdaG');
errorL2 = rmfield(errorL2,'lambda');
errorL2 = rmfield(errorL2,'SN05');
%
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
symbols = ['*'; '+'; 'o'; 'x'; 'd'];
AT = newTimeStepSizeVector;
slopes = zeros(numberOfVariables, 1);
for markerCounter = 1:numberOfVariables
    % Ausgleichsgeraden
    % s = scatter(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))));
    % s.Marker = 'x';
    yp = polyfit(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))), 1);
    % plot(log10(AT), polyval(yp, log10(AT)));
    slopes(markerCounter) = yp(1);
    plot(log10(AT), log10(vertcat(errorL2.(variableNames{markerCounter}))), strcat('-', symbols(markerCounter)));
    % loglog(AT,vertcat(errorL2.(variableNames{markerCounter})),strcat('-',symbols(markerCounter)));
end
axis tight;
grid on;
% set(gca,'xscale','log','yscale','log')
title('Log_{10} of L2-Norm of error');
xlabel('\Delta t');
ylabel('error');
% legendEntries = cell(numberOfVariables, 1);
legendEntries = variableNames;
% for ii=1:numberOfVariables
%     switch variableNames{ii}
%         case 'phi'
%             legendEntries{ii} = '\varphi';
%         case 'C'
%             legendEntries{ii} = '\vec{C}';
%         case 'G'
%             legendEntries{ii} = '\vec{G}';
%         case 'J'
%             legendEntries{ii} = 'J';
%         case 'lambdaC'
%             legendEntries{ii} = '\vec{\Lambda}^{\subC}';
%         case 'lambdaG'
%             legendEntries{ii} = '\vec{\Lambda}^{\subG}';
%         case 'lambda'
%             legendEntries{ii} = '\lambda';
%         case 'SN1'
%             legendEntries{ii} = '\vec{S}_{n+1}';
%         case 'SN05'
%             legendEntries{ii} = '\vec{S}_{n+1/2}';
%     end
% end
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

matlab2tikz([fileName 'timeConvergenceDisplemcent','.tikz'],'width','\figW','height','\figH')

function s = autoNum2Str(x)
    s = sprintf('%.15f', x);           % genügend Nachkommastellen
    s = regexprep(s, '0+$', '');       % entferne trailing zeros
    s = regexprep(s, '\.$', '');       % entferne Punkt am Ende
end

