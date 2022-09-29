import numpy as np
from joblib import Parallel, delayed
from commonFunctions.expandEdof import expandEdof
from commonFunctions.computeDeterminant import computeDeterminant
#from continuumClasses.solidClass.elementLoop import elementLoop

def displacementParallelHookeEndpoint(obj,setupObject):
    # Creates the residual and the tangent of the given obj.
    #
    # Syntax
    #
    # out = linearEndPoint(obj,'PropertyName',PropertyValue)
    #
    # Description
    #
    # homogenous linear-elastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
    #
    # 08.02.2015 M.FRANKE
    
    # Check input
    
    # Shape functions
    globalFullEdof = obj.globalFullEdof
    edof = obj.edof
    numberOfGausspoints = obj.numberOfGausspoints
    gaussWeight = obj.shapeFunctions['gaussWeight']
    NAll = obj.shapeFunctions['N'].T
    dNrAll = obj.shapeFunctions['dNr'].T
    qR = obj.qR
    qN1 = obj.qN1
    numberOfElements = globalFullEdof.shape[0]
    numberOfDOFs = globalFullEdof.shape[1]
    dimension = obj.dimension
    
    Lambda = obj.materialData['Lambda']
    mu = obj.materialData['mu']
    DMat = np.array([[Lambda+2*mu, Lambda, Lambda, 0, 0, 0],
                     [Lambda, Lambda+2*mu, Lambda, 0, 0, 0],
                     [Lambda, Lambda, Lambda+2*mu, 0, 0, 0],
                     [0, 0, 0, mu, 0, 0],
                     [0, 0, 0, 0, mu, 0],
                     [0, 0, 0, 0, 0, mu]])
    
    # Create residual and tangent
    dataFE = {'edofE': [None]*numberOfElements,
              'Re': [None]*numberOfElements,
              'ePot': [None]*numberOfElements,
              'pI': [None]*numberOfElements,
              'pJ':[None]*numberOfElements,
              'pK':[None]*numberOfElements}
   
    # dataParallel = Parallel(n_jobs=2)(delayed(elementLoop)(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat, dataFE) for e in np.arange(0,numberOfElements))
    dataParallel = Parallel(n_jobs=24)(delayed(element)(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat) for e in range(numberOfElements))
    for e in np.arange(0, np.size(dataParallel)):
        dataFE['Re'][e] = dataParallel[e]['Re']
        dataFE['edofE'][e] = dataParallel[e]['edofE']        
        dataFE['ePot'][e] = dataParallel[e]['ePot']        
        dataFE['pI'][e] = dataParallel[e]['pI']
        dataFE['pJ'][e] = dataParallel[e]['pJ']
        dataFE['pK'][e] = dataParallel[e]['pK']        
        
    return dataFE
        
# def elementLoop(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat):        
#     dataFE = element(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat) 
#     return dataFE

def element(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat):
    dataFE = {}
    edofH1,edofH2 = expandEdof(globalFullEdof[e,:])
    dataFE['pI'] = edofH1.T
    dataFE['pJ'] = edofH2.T
    dataFE['edofE'] = globalFullEdof[e,:].T  # evtl. noch konvertieren in Double notwendig?!
    # Element routine
    Re = np.zeros(numberOfDOFs)#.reshape(-1,1)
    Ke = np.zeros((numberOfDOFs, numberOfDOFs))
    ePot = 0
    edN1 = qN1[edof[e,:],].T
    edRef = qR[edof[e,:],].T
    uN1 = edN1 - edRef
    J = qR[edof[e,:],].T @ dNrAll
    # Run through all Gauss points
    for k in range(numberOfGausspoints):
        indx = range(dimension * k, dimension*(k+1))
        detJ = computeDeterminant(J[:,indx].T)
        if (detJ < 10*np.spacing(1)):
            print('Jacobi determinant equal or less than zero.'),  exit()

        dNx = np.linalg.solve(J[:,indx].T,dNrAll[:,indx].T)
        B = np.zeros((6,numberOfDOFs))
        B[0,0::3] = dNx[0,:]
        B[1,1::3] = dNx[1,:]
        B[2,2::3] = dNx[2,:]
        B[3,0::3] = dNx[1,:]
        B[3,1::3] = dNx[0,:]
        B[4,1::3] = dNx[2,:]
        B[4,2::3] = dNx[1,:]
        B[5,0::3] = dNx[2,:]
        B[5,2::3] = dNx[0,:]
        Re = Re + uN1.T.reshape(-1) @ (B.T @ DMat @ B) * detJ * gaussWeight[k]
        Ke = Ke + B.T @ DMat @ B * detJ * gaussWeight[k]

    dataFE['Re'] = Re
    dataFE['pK'] = Ke.T.reshape(-1)
    dataFE['ePot'] = ePot
    return dataFE