function dofObject = runNewton(setupObject,dofObject)
%% initialization
dofObject.initializeGlobalDofs('initializeGlobalDofs');
dofObject.initializeQV('initializeQV');
dofObject.initializeShapeFunctions('initializeShapeFunctions');
% global states
initialize(dofObject,{'qN','vN'});
updateGlobalField(dofObject,{'qN','vN'});
% compute the constant mass matrix
M = dofObject.callMassMatrixElement;
% compute solveDof from globalDirichletList
solveDof = dofObject.solveDof;
% initial guess
dofObject.qN1 = dofObject.qN;
dofObject.vN1 = dofObject.vN;

%% time loop
time = 0;
DT = setupObject.timeStepSize;
for timeStep = 1:setupObject.totalTimeSteps
    fprintf('|---------- time step %4.5d of %4.5d ----------| \n',timeStep,setupObject.totalTimeSteps);
    time = time + DT;
    % update fields & objects for every timestep
    dofObject.qN = dofObject.qN1;
    dofObject.vN = dofObject.vN1;
    updateTimeDependentField(dofObject,{'qN1'},time,'updateTimeDependentField');
    updateContinuumField(dofObject,{'qN','vN','qN1'},'updateContinuumFieldPreNewtonLoop')
    
    %% Newton loop
    newtonStep = 0;
    while newtonStep <= setupObject.newtonMaximumSteps
        newtonStep = newtonStep + 1;
        dataFE = callElements(dofObject,setupObject,'callElements');
        R = sparse(vertcat(dataFE(:).edofE),1,vertcat(dataFE(:).Re),dofObject.totalNumberOfDofs,1);
        R = setupObject.factorIntegrator(1)/DT^2*M*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator(2)/DT*M*dofObject.vN + R(1:dofObject.totalNumberOfDofs);
        normR = norm(R(solveDof));
        fprintf('|      Iteration %2d: norm(R) = %9.6e    |\n',newtonStep,normR)
        if (newtonStep > setupObject.newtonMaximumSteps || normR > setupObject.newtonMaximumResiduumNorm)
            error('Newton iteration failed!');
        elseif normR < setupObject.newtonTolerance
            % converged state
            newtonStep = setupObject.newtonMaximumSteps + 1;
        else
            K = sparse(vertcat(dataFE(:).pI), vertcat(dataFE(:).pJ), vertcat(dataFE(:).pK),dofObject.totalNumberOfDofs,dofObject.totalNumberOfDofs);
            K = K + setupObject.factorIntegrator(1)/DT^2*M;
            % direct solver
            dofObject.qN1(solveDof,1) = dofObject.qN1(solveDof,1) - K(solveDof,solveDof)\R(solveDof,1);
            % update objects for every Newton iteration
            updateContinuumField(dofObject,{'qN1'},'updateContinuumFieldNewtonLoop')
        end
    end
    
    dofObject.vN1 = setupObject.factorIntegrator(1)/DT*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator(2)*dofObject.vN;
    plotScript
end
end