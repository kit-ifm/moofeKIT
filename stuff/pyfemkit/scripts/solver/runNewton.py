import numpy as np
from sys import exit
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
#from scikits.umfpack import spsolve
from commonFunctions.callElements import callElements

def runNewton(setupObject,dofObject):
    # initialization
    dofObject.initializeGlobalDofs('initializeGlobalDofs')
    dofObject.initializeQV('initializeQV')
    dofObject.initializeShapeFunctions('initializeShapeFunctions')
    # global states
    dofObject.initialize(["qN","vN"]);
    dofObject.updateGlobalField( ["qN","vN"]);
    # compute the constant mass matrix
    M = dofObject.callMassMatrixElement
    # compute solveDof from globalDirichletList
    solveDof = dofObject.solveDof()
    # initial guess
    dofObject.qN1 = dofObject.qN.copy()
    dofObject.vN1 = dofObject.vN.copy()
    # time loop
    time = 0
    DT = setupObject.timeStepSize
    for timeStep in range(1, setupObject.totalTimeSteps + 1):
        print('|---------- time step %4.5d of %4.5d ----------|' % (timeStep,setupObject.totalTimeSteps))
        time += DT
        # update fields & objects for every timestep
        dofObject.qN = dofObject.qN1.copy()
        dofObject.vN = dofObject.vN1.copy()
        dofObject.updateTimeDependentField(["qN1"],time,'updateTimeDependentField')
        dofObject.updateContinuumField(["qN","vN","qN1"],'updateContinuumFieldPreNewtonLoop')
    
        # Newton loop
        newtonStep = 0
        while newtonStep <= setupObject.newtonMaximumSteps:
            newtonStep = newtonStep + 1
            dataFE = callElements(dofObject,setupObject,'callElements')
            R = csr_matrix((np.concatenate(dataFE['Re']),(np.concatenate(dataFE['edofE']), np.zeros(np.size(np.concatenate(dataFE['edofE'])), dtype="i"))))
            R = csr_matrix(setupObject.factorIntegrator[0]/DT**2*M*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator[1]/DT*M*dofObject.vN) + R.reshape(1,-1)
            normR = np.linalg.norm(R[0,solveDof].toarray())
            print('|      iteration %2d: norm(R) = %9.6e    |' % (newtonStep,normR))
            if (newtonStep > setupObject.newtonMaximumSteps or normR > setupObject.newtonMaximumResiduumNorm):
                exit('Newton iteration failed!')
            elif (normR < setupObject.newtonTolerance):
                # converged state
                newtonStep = setupObject.newtonMaximumSteps + 1;
            else:
                K = csr_matrix((np.concatenate(dataFE['pK']),(np.concatenate(dataFE['pI']), np.concatenate(dataFE['pJ'])))) 
                K = K + setupObject.factorIntegrator[0]/DT**2*M
                # direct solver
                dofObject.qN1[solveDof] = dofObject.qN1[solveDof] - spsolve(K[np.ix_(solveDof,solveDof)], R[0,solveDof].T)
                # update objects for every Newton iteration
                dofObject.updateContinuumField(["qN1"],'updateContinuumFieldNewtonLoop')
        dofObject.vN1 = setupObject.factorIntegrator[0]/DT*(dofObject.qN1 - dofObject.qN) - setupObject.factorIntegrator[1]*dofObject.vN
        # plotScript
            
    return dofObject            
       
