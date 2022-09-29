function [rData, kData, elementEnergy, array] = deadLoad1D(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% load objects
storageFEObject = obj.storageFEObject;
globalEdof = obj.globalEdof;
% tol = 1e-12; % numerische Toleranz für if-Bedingung
if (setupObject.integrator == 'Midpoint')
    % für Mittelpunktsregel: Auswertung bei t_(N + 1/2)
    DT = setupObject.timeStepSize;
    forceVector = obj.forceVector*obj.timeFunction(obj.time - 0.5*DT);
else
    forceVector = obj.forceVector*obj.timeFunction(obj.time);
end
timeStep = setupObject.timeStep;
elementEnergy.externalEnergy = 0;
rData{1} = rData{1} - forceVector;
obj.nodalForce(size(globalEdof,1)) = forceVector;
% elementDataFE = storageFEObject.assignElementDataFE(initialDataFE, array, globalEdof, e, dimension, [], [], computePostData);
% dataFE(e) = elementDataFE;
% storageFEObject.storeDataFE(dataFE, computePostData);
% FIXME: include externalEnergy
% obj.ePot(setupObject.timeStep).externalEnergy = externalEnergy;
end