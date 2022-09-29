if 1
    if ~isempty(dofObject.dirichletDof)
        doSolve = true(dofObject.totalNumberOfDofs,1);
        doSolve(dofObject.dirichletDof) = false;
        delta = zeros(dofObject.totalNumberOfDofs,1);
        % FIXME
        if newtonStep == 1 
            delta(doSolve) = K(doSolve,doSolve)\(-dofObject.R(doSolve) - K(doSolve,~doSolve)*(dofObject.qN1(~doSolve)-dofObject.qN(~doSolve)));
        else
            delta(doSolve) = K(doSolve,doSolve)\(-dofObject.R(doSolve));
        end
    else
        delta = -K(solveDof,solveDof)\dofObject.R(solveDof,1);
    end
    dofObject.qN1(doSolve) = dofObject.qN1(doSolve) + delta(doSolve);
else
    delta = -K(solveDof,solveDof)\dofObject.R(solveDof,1);
    dofObject.qN1(solveDof,1) = dofObject.qN1(solveDof,1) + delta;
end

