%% Script for plotting the total energy
clear;
clc;

momentum = 'totalLinearMomentum'; % linearMomentum, totalLinearMomentum, angularMomentum, totalAngularMomentum
filename = 'rotatingXElectroThermoMechanics';
lastTimeStep = 150;
exportAsTikz = true;

%% plot momentum
load(strcat(filename, num2str(lastTimeStep), '.mat'));

timeVector = getTime(dofObject.postDataObject,setupObject);
if strcmpi(momentum, 'linearMomentum') || strcmpi(momentum, 'totalLinearMomentum')
    [linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
    if strcmpi(momentum, 'linearMomentum')
        momentumData = linearMomentum;
    elseif strcmpi(momentum, 'totalLinearMomentum')
        momentumData = totalLinearMomentum;
    end
elseif strcmpi(momentum, 'angularMomentum') || strcmpi(momentum, 'totalAngularMomentum')
    [angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
    if strcmpi(momentum, 'angularMomentum')
        momentumData = angularMomentum;
    elseif strcmpi(momentum, 'totalAngularMomentum')
        momentumData = totalAngularMomentum;
    end
end

figure;
plot(timeVector, momentumData);

%% export total energy as tikz
if exportAsTikz
    matlab2tikz(strcat(filename, upper(momentum(1)), momentum(2:end), '.tikz'), 'width', '\fwidth', 'height', '\fheight');
end