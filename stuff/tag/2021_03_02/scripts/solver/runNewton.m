% TODO: entweder automatisieren oder evtl. erstmal weglassen, da man sich
% auf Konvention einigen kann: solidObject(1) = solidClass; solidObject(2)
% = solidClass; etc.

% initialize Setup (assign global dofs, initialize states q,v={ref,n,n+1}, 
totalNumberOfDofs = 0;
totalNumberOfDofs = assignGlobalDofs(solidObject,totalNumberOfDofs); 
initializeQV(solidObject); 
assignShapefunctions(solidObject)
% TODO: datamanager

%% constant mass matrices
massData = massMatrix(solidObject);
% sparse assembly
M = sparse(vertcat(massData(:).indexMi),vertcat(massData(:).indexMj),vertcat(massData(:).MeVector),totalNumberOfDofs,totalNumberOfDofs);
clear massData

%% initialize global states
qN = zeros(totalNumberOfDofs,1);
vN = zeros(totalNumberOfDofs,1);
% TODO: automatisieren fuer alle systeme
for index = 1:numel(solidObject)
    qN(solidObject(index).globalNodesDOF) = solidObject(index).QN; 
    vN(solidObject(index).globalNodesDOF) = solidObject(index).VN;
end
% initial guess
qN1 = qN;

dirichletDofs = unique(vertcat(dirichletObject.masterGlobalNodesDof));
solveDof = 1:totalNumberOfDofs;
solveDof(dirichletDofs(:,1)) = [];

% integrator settings
if strcmpi(setupObject.integrator,'endpoint')
    factorIntegrator = 1;
    factorIntegratorVelocity = 0;
elseif (strcmpi(setupObject,'midpoint') || strcmpi(setupObject,'discreteGradient'))
    factorIntegrator = 2;
    factorIntegratorVelocity = 0;
end

%% time loop
time = 0;
DT = setupObject.timeStepSize;
for timeStep = 1:setupObject.totalTimeSteps
    fprintf('|---------- time step %4.5d of %4.5d ----------| \n',timeStep,setupObject.totalTimeSteps);
    time = time + DT;
    % update fields & objects
    qN = qN1;
    updateField(solidObject,'QN',qN);
    updateField(solidObject,'VN',vN);
    % update time dependent objects
    for index = 1:numel(dirichletObject)
        dirichletObject(index).time = time; %#ok<SAGROW>
        qN1(dirichletObject(index).masterGlobalNodesDof) = dirichletObject(index).qN1;
    end
    updateField(solidObject,'QN1',qN1);        
    % Newton loop
    newtonStep = 0;
    while newtonStep <= setupObject.newtonMaximumSteps
        newtonStep = newtonStep + 1;
        % TODO: for all continuum Objects
        dataFE = callElements(solidObject,setupObject);
        % sparse assembling residuum
        R = sparse(vertcat(dataFE(:).edofE),1,vertcat(dataFE(:).Re),totalNumberOfDofs,1);
        R = factorIntegrator/DT^2*M*(qN1 - qN) - factorIntegrator/DT*M*vN + R(1:totalNumberOfDofs);
        normR = norm(R(solveDof));
        fprintf('Iteration %2d: norm(R) = %9.6e \n',newtonStep,normR)
        if (newtonStep > setupObject.newtonMaximumSteps || normR > setupObject.newtonMaximumResiduumNorm)
            error('Newton iteration failed!');
        elseif normR < setupObject.newtonTolerance
            % converged state
            newtonStep = setupObject.newtonMaximumSteps + 1;
        else
            % sparse assembling tangent
            K = sparse(vertcat(dataFE(:).pI), vertcat(dataFE(:).pJ), vertcat(dataFE(:).pK),totalNumberOfDofs,totalNumberOfDofs);
            K = K + factorIntegrator/DT^2*M;
            % direct solver
            qN1(solveDof,1) = qN1(solveDof,1) - K(solveDof,solveDof)\R(solveDof,1);            
            % update
            updateField(solidObject,'QN1',qN1);
        end
    end
        
    % update field & objects
    vN = factorIntegrator/DT*(qN1 - qN) - factorIntegratorVelocity*vN;
    updateField(solidObject,'VN',vN);
    for index = 1:numel(solidObject)
        solidObject(index).QN = solidObject(index).QN1;
    end
    
%     plotScript
end
