%% initialization
dofsObject.collectSolverObjects;
dofsObject.initializeGlobalDofs('initializeGlobalDofs');
dofsObject.initializeQV('initializeQV'); 
dofsObject.initializeShapeFunctions('initializeShapeFunctions');
% global states
initialize(dofsObject,{'qN','vN'});
updateGlobalField(dofsObject,{'qN','vN'});
% compute the constant mass matrix
M = dofsObject.callMassMatrixElement;
% compute solveDof from globalDirichletList
solveDof = dofsObject.solveDof;
% initial guess
dofsObject.qN1 = dofsObject.qN;
dofsObject.vN1 = dofsObject.vN;

%% time loop
time = 0;
DT = setupObject.timeStepSize;
for timeStep = 1:setupObject.totalTimeSteps
    fprintf('|---------- time step %4.5d of %4.5d ----------| \n',timeStep,setupObject.totalTimeSteps);
    time = time + DT;
    % update fields & objects for every timestep
    dofsObject.qN = dofsObject.qN1; 
    dofsObject.vN = dofsObject.vN1;
    updateTimeDependentField(dofsObject,{'qN1'},time,'updateTimeDependentField');
    updateContinuumField(dofsObject,{'qN','vN','qN1'},'updateContinuumFieldPreNewtonLoop')

    %% Newton loop
    newtonStep = 0;
    while newtonStep <= setupObject.newtonMaximumSteps
        newtonStep = newtonStep + 1;
        dataFE = callElements(dofsObject,setupObject,'callElements');
        R = sparse(vertcat(dataFE(:).edofE),1,vertcat(dataFE(:).Re),dofsObject.totalNumberOfDofs,1);
        R = setupObject.factorIntegrator(1)/DT^2*M*(dofsObject.qN1 - dofsObject.qN) - setupObject.factorIntegrator(2)/DT*M*dofsObject.vN + R(1:dofsObject.totalNumberOfDofs);
        normR = norm(R(solveDof));
        fprintf('|      Iteration %2d: norm(R) = %9.6e    |\n',newtonStep,normR)
        if (newtonStep > setupObject.newtonMaximumSteps || normR > setupObject.newtonMaximumResiduumNorm)
            error('Newton iteration failed!');
        elseif normR < setupObject.newtonTolerance
            % converged state
            newtonStep = setupObject.newtonMaximumSteps + 1;
        else
            K = sparse(vertcat(dataFE(:).pI), vertcat(dataFE(:).pJ), vertcat(dataFE(:).pK),dofsObject.totalNumberOfDofs,dofsObject.totalNumberOfDofs);
            K = K + setupObject.factorIntegrator(1)/DT^2*M;
            % direct solver
            dofsObject.qN1(solveDof,1) = dofsObject.qN1(solveDof,1) - K(solveDof,solveDof)\R(solveDof,1);            
            % update objects for every Newton iteration
            updateContinuumField(dofsObject,{'qN1'},'updateContinuumFieldNewtonLoop')
        end
    end
    
    dofsObject.vN1 = setupObject.factorIntegrator(1)/DT*(dofsObject.qN1 - dofsObject.qN) - setupObject.factorIntegrator(2)*dofsObject.vN;    
    plotScript
end