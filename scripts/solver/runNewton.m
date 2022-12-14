function dofObject = runNewton(setupObject,dofObject)
%% initialization
% initialize saveObject
initialize(setupObject.saveObject, setupObject);
% initialize continuumObjects
dofObject.initializeObjectsShapeFunctions;
dofObject.initializeObjectsQV;
dofObject.initializeObjectsGlobalDofs;
% global states
initialize(dofObject,{'qR','qN','vN'});
updateGlobalField(dofObject,{'qR','qN','vN'});
% compute the constant mass matrix
M = dofObject.callMassMatrixElement(setupObject);
% initialize postDataObject
initialize(dofObject.postDataObject,dofObject,M);
% compute solveDof from globalDirichletList
doSolve = dofObject.doSolve;
% initial guess
dofObject.qN1 = dofObject.qN;
dofObject.vN1 = dofObject.vN;
delta = zeros(size(dofObject.qN,1),1);
%% time loop
DT = setupObject.timeStepSize;
for timeStep = 1:setupObject.totalTimeSteps
    setupObject.timeStep = timeStep;
    fprintf('|--- time step %4.5d of %4.5d; time %4.4f ---| \n',timeStep,setupObject.totalTimeSteps,setupObject.time);
    setupObject.time = timeStep*DT;
    % update fields & objects for every timestep
    dofObject.qN = dofObject.qN1;
    dofObject.vN = dofObject.vN1;
    dofObject.updateObjectsContinuumFieldPreNewtonLoop(setupObject,{'qN','vN','qN1'})
    updateTimeDependentFieldPreNewtonLoop(dofObject,'qN1',setupObject.time);
    updateHistoryField(dofObject);
    %% Newton loop
    setupObject.newton.step = 0;
    while setupObject.newton.step <= setupObject.newton.maximumSteps
        setupObject.newton.step = setupObject.newton.step + 1;
        dataFE = [];
        for index1 = 1:dofObject.numberOfContinuumObjects
            continuumObject = dofObject.listContinuumObjects{index1};
            if continuumObject.callElements
                dataFEContinuumObject = callElements(continuumObject,setupObject,false);
                dataFE = [dataFE, dataFEContinuumObject];
            end
        end
        % assembly via sparse command, for more details see https://arxiv.org/pdf/1305.3122.pdf
        dofObject.R = sparse(vertcat(dataFE(:).indexReI),1,vertcat(dataFE(:).Re),dofObject.totalNumberOfDofs,1);
        dofObject.R = M*(setupObject.factorIntegrator(1)/DT^2*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator(2)/DT*dofObject.vN) + dofObject.R(1:dofObject.totalNumberOfDofs);
        K = sparse(vertcat(dataFE(:).indexKeI), vertcat(dataFE(:).indexKeJ), vertcat(dataFE(:).Ke),dofObject.totalNumberOfDofs,dofObject.totalNumberOfDofs);
        K = K + setupObject.factorIntegrator(1)/DT^2*M;
        if setupObject.newton.step == 1
            % TODO @Marlon: documentation!
            updateTimeDependentFieldNewtonLoop(dofObject,'qN1',setupObject.time);
            dofObject.R(doSolve) = dofObject.R(doSolve) + K(doSolve,~doSolve)*(dofObject.qN1(~doSolve)-dofObject.qN(~doSolve));
            % TODO @Marlon: documentation!
        end
        normR = norm(dofObject.R(doSolve));
        if setupObject.newton.enforceIteration && (setupObject.newton.step == 1) && (normR < setupObject.newton.tolerance)
            fprintf('|      iteration %2d: norm(R) = %9.5e    |\n',0,normR)
            normR = 1e5;
        end
        fprintf('|      iteration %2d: norm(R) = %9.5e    |\n',setupObject.newton.step,normR)
        if (setupObject.newton.step > setupObject.newton.maximumSteps || normR > setupObject.newton.maximumResiduumNorm)
            error('Newton iteration failed!');
        elseif normR < setupObject.newton.tolerance
            % converged state
            break;
        else
            % direct solver
            delta(doSolve) = -solverLinearSystem(setupObject, K(doSolve,doSolve), dofObject.R(doSolve));
            dofObject.qN1(doSolve) = dofObject.qN1(doSolve) + delta(doSolve);
            % update objects for every Newton iteration
            updateObjectsContinuumFieldNewtonLoop(dofObject,setupObject,delta,'qN1')
        end
    end
    dofObject.vN1 = setupObject.factorIntegrator(1)/DT*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator(3)*dofObject.vN;
    updateObjectsContinuumFieldPostNewtonLoop(dofObject,setupObject,'vN1');
    % store momentum maps and kinetic energy
    update(dofObject.postDataObject,dofObject,timeStep,setupObject.time,M); 
    % plot fe body
    plotScript
    % save data as .mat file
        saveDataAsMatFile(setupObject.saveObject, setupObject, dofObject, timeStep);
    end

% terminate saveObject
terminate(setupObject.saveObject);
end