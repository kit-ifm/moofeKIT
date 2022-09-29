import numpy as np
# from scipy import linalg, kron
from commonFunctions.expandEdof import expandEdof
from commonFunctions.kron_A_dimension import kron_A_dimension
from commonFunctions.computeDeterminant import computeDeterminant

def massMatrixElement(self):
    globalFullEdof = self.globalFullEdof
    edof = self.edof
    numberOfGausspoints = self.numberOfGausspoints
    gaussWeight = self.shapeFunctions['gaussWeight']
    NAll = self.shapeFunctions['N'].T
    dNrAll = self.shapeFunctions['dNr'].T
    rho = self.materialData['rho']
    qR = self.qR
    numberOfElements = globalFullEdof.shape[0]
    dimension = self.dimension
    # I = np.eye(dimension)
    # additionalFields = self.additionalFields
    
    massDataFE = {'MeVector':[None]*numberOfElements,
              'indexMi':[None]*numberOfElements,
              'indexMj':[None]*numberOfElements}
    
    # TODO: Parallelisierung
    for e in np.arange(0,numberOfElements):
        globalEdofH1, globalEdofH2 = expandEdof(globalFullEdof[e,:])
        massDataFE['indexMi'][e] = globalEdofH1
        massDataFE['indexMj'][e] = globalEdofH2
        numberOfDofs = len(globalFullEdof[e,:])
        Me = np.zeros((numberOfDofs, numberOfDofs))
        J = qR[edof[e,:],].T @ dNrAll
        for gp in range(numberOfGausspoints):
            index2 = range(dimension * gp, dimension*(gp+1))
            # detJ = linalg.det(J[:,index2].T,True,False)
            detJ = computeDeterminant(J[:,index2].T)
            if detJ <=0:
                print('Jacobideterminant less or equal zero.')
            N = NAll[gp,:]
            # Nkron = kron(N,I)
            Nkron = kron_A_dimension(N.reshape(1,-1),dimension)
            Me = Me + rho * Nkron.T @ Nkron * detJ * gaussWeight[gp]
        massDataFE['MeVector'][e] = Me.reshape(np.size(Me))

    return massDataFE