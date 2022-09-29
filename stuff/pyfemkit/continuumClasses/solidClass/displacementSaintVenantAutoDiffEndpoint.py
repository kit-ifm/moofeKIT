import auto_diff
import numpy as np
from commonFunctions.expandEdof import expandEdof
from commonFunctions.computeDeterminant import computeDeterminant

def displacementSaintVenantAutoDiffEndpoint(obj,setupObject):
    # Creates the residual and the tangent of the given obj.
    #
    # Syntax
    #
    # out = stVenantEndPoint(obj,'PropertyName',PropertyValue)
    #
    # Description
    #
    # St. Venant Kirchhoff strain-energy function, evaluated at time n+1, i.e. implicid euler method.
    #
    # 10.01.2012 C.HESCH
    
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
    
    # Create residual and tangent
    Lambda = obj.materialData['Lambda']
    mu = obj.materialData['mu']
    DMat = np.array([[Lambda+2*mu, Lambda, Lambda, 0, 0, 0],
                     [Lambda, Lambda+2*mu, Lambda, 0, 0, 0],
                     [Lambda, Lambda, Lambda+2*mu, 0, 0, 0],
                     [0, 0, 0, mu, 0, 0],
                     [0, 0, 0, 0, mu, 0],
                     [0, 0, 0, 0, 0, mu]])
    
    I = np.eye(dimension)
    # Create residual and tangent
    out = {'edofE': [None]*numberOfElements,
              'Re':[None]*numberOfElements,
              'ePot':[None]*numberOfElements,
              'pI': [None]*numberOfElements,
              'pJ':[None]*numberOfElements,
              'pK':[None]*numberOfElements}
    
    # TODO: Parallelisierung
    for e in np.arange(0,numberOfElements):
        edofH1,edofH2 = expandEdof(globalFullEdof[e,:])
        out['pI'][e] = edofH1.T
        out['pJ'][e] = edofH2.T
        out['edofE'][e] = globalFullEdof[e,:].T  # evtl. noch konvertieren in Double notwendig?!
        # Element routine
        Re = np.zeros(numberOfDOFs)
        Ke = np.zeros((numberOfDOFs, numberOfDOFs))
        ePot = 0
        edN1 = qN1[edof[e,:],].T
        q = np.array([edN1.T.reshape(-1)]).T
        # J = qR[np._ix(edof[e,:],range(dimension))].T @ dNrAll
        J = qR[edof[e,:],].T @ dNrAll
        # Run through all Gauss points
        for k in range(numberOfGausspoints):
            indx = range(dimension * k, dimension*(k+1))
            # detJ = np.linalg.det(J[:,indx].T)
            detJ = computeDeterminant(J[:,indx].T)
            if (detJ < 10*np.spacing(1)):
                print('Jacobi determinant equal or less than zero.'),  exit()
            # Strain energy
            # ePot = ePot + 0.5*En1_v.T@DMat@En1_v*detJ*gaussWeight[k]
            # re,ke = residualFunction(edN1,J,indx,dNrAll,numberOfDOFs,DMat,I,dimension)
            re = residualFunction(q,J,indx,dNrAll,numberOfDOFs,DMat,I,dimension)
            with auto_diff.AutoDiff(q) as z:
                re_eval = residualFunction(z,J,indx,dNrAll,numberOfDOFs,DMat,I,dimension)
                ke = auto_diff.jacobian(re_eval)
            Re = Re + re.flatten() * detJ * gaussWeight[k]
            Ke = Ke + ke * detJ * gaussWeight[k]

        out['Re'][e] = Re
        out['pK'][e] = Ke.T.reshape(-1)
        out['ePot'][e] = ePot
        
    return out

def residualFunction(q,J,indx,dNrAll,numberOfDOFs,DMat,I,dimension):
    # edN1 = q.reshape(8,3).T
    edN1 = np.zeros((3,8))
    edN1[0,0] = q[0,0]
    edN1[1,0] = q[1,0]
    edN1[2,0] = q[2,0]
    edN1[0,1] = q[3,0]
    edN1[1,1] = q[4,0]
    edN1[2,1] = q[5,0]
    edN1[0,2] = q[6,0]
    edN1[1,2] = q[7,0]
    edN1[2,2] = q[8,0]
    edN1[0,3] = q[9,0]
    edN1[1,3] = q[10,0]
    edN1[2,3] = q[11,0]
    edN1[0,4] = q[12,0]
    edN1[1,4] = q[13,0]
    edN1[2,4] = q[14,0]
    edN1[0,5] = q[15,0]
    edN1[1,5] = q[16,0]
    edN1[2,5] = q[17,0]
    edN1[0,6] = q[18,0]
    edN1[1,6] = q[19,0]
    edN1[2,6] = q[20,0]
    edN1[0,7] = q[21,0]
    edN1[1,7] = q[22,0]
    edN1[2,7] = q[23,0]    
    dNx = np.linalg.solve(J[:,indx].T,dNrAll[:,indx].T)
    FAkt = edN1 @ dNx.T
    # B-matrix
    BAkt = np.zeros((6,numberOfDOFs))
    BAkt[0,0::3]=FAkt[0,0]*dNx[0,:]
    BAkt[0,1::3]=FAkt[1,0]*dNx[0,:]
    BAkt[0,2::3]=FAkt[2,0]*dNx[0,:]
    BAkt[1,0::3]=FAkt[0,1]*dNx[1,:]
    BAkt[1,1::3]=FAkt[1,1]*dNx[1,:]
    BAkt[1,2::3]=FAkt[2,1]*dNx[1,:]
    BAkt[2,0::3]=FAkt[0,2]*dNx[2,:]
    BAkt[2,1::3]=FAkt[1,2]*dNx[2,:]
    BAkt[2,2::3]=FAkt[2,2]*dNx[2,:]
    BAkt[3,0::3]=FAkt[0,0]*dNx[1,:]+FAkt[0,1]*dNx[0,:]
    BAkt[3,1::3]=FAkt[1,0]*dNx[1,:]+FAkt[1,1]*dNx[0,:]
    BAkt[3,2::3]=FAkt[2,0]*dNx[1,:]+FAkt[2,1]*dNx[0,:]
    BAkt[4,0::3]=FAkt[0,1]*dNx[2,:]+FAkt[0,2]*dNx[1,:]
    BAkt[4,1::3]=FAkt[1,1]*dNx[2,:]+FAkt[1,2]*dNx[1,:]
    BAkt[4,2::3]=FAkt[2,1]*dNx[2,:]+FAkt[2,2]*dNx[1,:]
    BAkt[5,0::3]=FAkt[0,0]*dNx[2,:]+FAkt[0,2]*dNx[0,:]
    BAkt[5,1::3]=FAkt[1,0]*dNx[2,:]+FAkt[1,2]*dNx[0,:]
    BAkt[5,2::3]=FAkt[2,0]*dNx[2,:]+FAkt[2,2]*dNx[0,:]
    # Cauchy-Green tensor
    CAkt = FAkt.T@FAkt
    # Green-Lagrange tensor
    En1   = 0.5*(CAkt-I)
    En1_v = np.array([[En1[0,0]], [En1[1,1]], [En1[2,2]], [2*En1[0,1]], [2*En1[2,1]], [2*En1[2,0]]])
    # Stresses
    DW1_v = 1/2*DMat@En1_v
    # DW1 = np.array([[DW1_v[0], DW1_v[3], DW1_v[5]], [DW1_v[3], DW1_v[1], DW1_v[4]], [DW1_v[5], DW1_v[4], DW1_v[2]]])
    # Residual
    re = 2*BAkt.T@DW1_v
    # Tangent
    # D2W1 = 1/4*DMat
    # A1 = 2*dNx.T@DW1@dNx
    # MAT = np.zeros((numberOfDOFs,numberOfDOFs))
    # for g in np.arange(0,dimension):
    #     MAT[np.ix_(range(g,numberOfDOFs,dimension),range(g,numberOfDOFs,dimension))] = A1
    # ke = 4*BAkt.T@D2W1@BAkt + MAT
    return re