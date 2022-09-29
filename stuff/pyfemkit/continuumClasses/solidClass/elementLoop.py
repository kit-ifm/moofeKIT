# from continuumClasses.solidClass.element import element
import numpy as np
# from scipy import linalg
from commonFunctions.expandEdof import expandEdof
from commonFunctions.computeDeterminant import computeDeterminant

# def elementLoop(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat):        
#     dataFE = element(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat) 
#     return dataFE

# def element(e, globalFullEdof, edof, numberOfGausspoints, gaussWeight, NAll, dNrAll, qR, qN1, numberOfElements, numberOfDOFs, dimension, DMat):
#     dataFE = {}
#     edofH1,edofH2 = expandEdof(globalFullEdof[e,:])
#     dataFE['pI'] = edofH1.T
#     dataFE['pJ'] = edofH2.T
#     dataFE['edofE'] = globalFullEdof[e,:].T  # evtl. noch konvertieren in Double notwendig?!
#     # Element routine
#     Re = np.zeros(numberOfDOFs)#.reshape(-1,1)
#     Ke = np.zeros((numberOfDOFs, numberOfDOFs))
#     ePot = 0
#     edN1 = qN1[edof[e,:],].T
#     edRef = qR[edof[e,:],].T
#     uN1 = edN1 - edRef
#     J = qR[edof[e,:],].T @ dNrAll
#     # Run through all Gauss points
#     for k in range(numberOfGausspoints):
#         indx = range(dimension * k, dimension*(k+1))
#         detJ = computeDeterminant(J[:,indx].T)
#         if (detJ < 10*np.spacing(1)):
#             print('Jacobi determinant equal or less than zero.'),  exit()

#         dNx = np.linalg.solve(J[:,indx].T,dNrAll[:,indx].T)
#         B = np.zeros((6,numberOfDOFs))
#         B[0,0::3] = dNx[0,:]
#         B[1,1::3] = dNx[1,:]
#         B[2,2::3] = dNx[2,:]
#         B[3,0::3] = dNx[1,:]
#         B[3,1::3] = dNx[0,:]
#         B[4,1::3] = dNx[2,:]
#         B[4,2::3] = dNx[1,:]
#         B[5,0::3] = dNx[2,:]
#         B[5,2::3] = dNx[0,:]
#         Re = Re + uN1.T.reshape(-1) @ (B.T @ DMat @ B) * detJ * gaussWeight[k]
#         Ke = Ke + B.T @ DMat @ B * detJ * gaussWeight[k]

#     dataFE['Re'] = Re
#     dataFE['pK'] = Ke.T.reshape(-1)
#     dataFE['ePot'] = ePot
#     return dataFE