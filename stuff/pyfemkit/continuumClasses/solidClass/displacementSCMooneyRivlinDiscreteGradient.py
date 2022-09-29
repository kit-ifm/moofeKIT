import numpy as np
from commonFunctions.expandEdof import expandEdof
from commonFunctions.computeDeterminant import computeDeterminant
from commonFunctions.wedgeA import wedgeA
from commonFunctions.wedgeAB import wedgeAB


def displacementSCMooneyRivlinDiscreteGradient(obj,setupObject):

    # %% Creates the residual and the tangent of the given obj.
    # %
    # % Syntax
    # %
    # % out = wedgeSC_mooneyRivlin_discreteGradient(obj,'PropertyName',PropertyValue)
    # %
    # % Description
    # %
    # % Mooney Rivlin strain-energy function, algorithmic stress formula, energs conserving.
    # %
    # % 20.09.2016 A.Janz
    
    globalFullEdof = obj.globalFullEdof
    edof = obj.edof
    numberOfGausspoints = obj.numberOfGausspoints
    gaussWeight = obj.shapeFunctions['gaussWeight']
    NAll = obj.shapeFunctions['N'].T
    dNrAll = obj.shapeFunctions['dNr'].T
    qR = obj.qR
    qN = obj.qN
    qN1 = obj.qN1
    numberOfElements = globalFullEdof.shape[0]
    numberOfDOFs = globalFullEdof.shape[1]
    dimension = obj.dimension
    
    # Create residual and tangent
    a = obj.materialData['a']
    b = obj.materialData['b']
    c = obj.materialData['c']
    d = obj.materialData['d']

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
        edN = qN[edof[e,:],].T
        edN05 = 0.5*(edN+edN1)
        # J = qR[np._ix(edof[e,:],range(dimension))].T @ dNrAll
        J = qR[edof[e,:],].T @ dNrAll
        # Run through all Gauss points
        for k in range(numberOfGausspoints):
            indx = range(dimension * k, dimension*(k+1))
            # detJ = np.linalg.det(J[:,indx].T)
            detJ = computeDeterminant(J[:,indx].T)
            if (detJ < 10*np.spacing(1)):
                print('Jacobi determinant equal or less than zero.'),  exit()
                
            dNx = np.linalg.solve(J[:,indx].T,dNrAll[:,indx].T)
            FN1 = edN1 @ dNx.T
            FN = edN @ dNx.T
            FN05 = edN05 @ dNx.T
            # B-matrix (midpoint configuration)
            BN05 = np.zeros((6,numberOfDOFs))
            BN05[0,0::3]=FN05[0,0]*dNx[0,:]
            BN05[0,1::3]=FN05[1,0]*dNx[0,:]
            BN05[0,2::3]=FN05[2,0]*dNx[0,:]
            BN05[1,0::3]=FN05[0,1]*dNx[1,:]
            BN05[1,1::3]=FN05[1,1]*dNx[1,:]
            BN05[1,2::3]=FN05[2,1]*dNx[1,:]
            BN05[2,0::3]=FN05[0,2]*dNx[2,:]
            BN05[2,1::3]=FN05[1,2]*dNx[2,:]
            BN05[2,2::3]=FN05[2,2]*dNx[2,:]
            BN05[3,0::3]=FN05[0,0]*dNx[1,:]+FN05[0,1]*dNx[0,:]
            BN05[3,1::3]=FN05[1,0]*dNx[1,:]+FN05[1,1]*dNx[0,:]
            BN05[3,2::3]=FN05[2,0]*dNx[1,:]+FN05[2,1]*dNx[0,:]
            BN05[4,0::3]=FN05[0,1]*dNx[2,:]+FN05[0,2]*dNx[1,:]
            BN05[4,1::3]=FN05[1,1]*dNx[2,:]+FN05[1,2]*dNx[1,:]
            BN05[4,2::3]=FN05[2,1]*dNx[2,:]+FN05[2,2]*dNx[1,:]
            BN05[5,0::3]=FN05[0,0]*dNx[2,:]+FN05[0,2]*dNx[0,:]
            BN05[5,1::3]=FN05[1,0]*dNx[2,:]+FN05[1,2]*dNx[0,:]
            BN05[5,2::3]=FN05[2,0]*dNx[2,:]+FN05[2,2]*dNx[0,:]
            BN05 = 2*BN05
            # B-matrix (current configuration)
            BN1 = np.zeros((6,numberOfDOFs))
            BN1[0,0::3]=FN1[0,0]*dNx[0,:]
            BN1[0,1::3]=FN1[1,0]*dNx[0,:]
            BN1[0,2::3]=FN1[2,0]*dNx[0,:]
            BN1[1,0::3]=FN1[0,1]*dNx[1,:]
            BN1[1,1::3]=FN1[1,1]*dNx[1,:]
            BN1[1,2::3]=FN1[2,1]*dNx[1,:]
            BN1[2,0::3]=FN1[0,2]*dNx[2,:]
            BN1[2,1::3]=FN1[1,2]*dNx[2,:]
            BN1[2,2::3]=FN1[2,2]*dNx[2,:]
            BN1[3,0::3]=FN1[0,0]*dNx[1,:]+FN1[0,1]*dNx[0,:]
            BN1[3,1::3]=FN1[1,0]*dNx[1,:]+FN1[1,1]*dNx[0,:]
            BN1[3,2::3]=FN1[2,0]*dNx[1,:]+FN1[2,1]*dNx[0,:]
            BN1[4,0::3]=FN1[0,1]*dNx[2,:]+FN1[0,2]*dNx[1,:]
            BN1[4,1::3]=FN1[1,1]*dNx[2,:]+FN1[1,2]*dNx[1,:]
            BN1[4,2::3]=FN1[2,1]*dNx[2,:]+FN1[2,2]*dNx[1,:]
            BN1[5,0::3]=FN1[0,0]*dNx[2,:]+FN1[0,2]*dNx[0,:]
            BN1[5,1::3]=FN1[1,0]*dNx[2,:]+FN1[1,2]*dNx[0,:]
            BN1[5,2::3]=FN1[2,0]*dNx[2,:]+FN1[2,2]*dNx[0,:]
            BN1 = 2*BN1
                            
            # Right Cauchy-Green tensor
            CN1 = FN1.T@FN1
            CN = FN.T@FN
            CN05 = 0.5*(CN+CN1);            #averaged
            
            # Cofactor
            GN1 = 0.5*wedgeA(CN1)
            GN = 0.5*wedgeA(CN)
            GN05 = 0.5*(GN+GN1)            #averaged
            GMid = 0.5*wedgeA(CN05)    #mid-configuration
            
            # Third invariant
            I3N1 = computeDeterminant(CN1)
            I3N = computeDeterminant(CN)
            I3N05 = 0.5*(I3N+I3N1)
            
            # Strain energy function
            ePot = ePot + (a*(np.trace(CN1)-3) + b*(np.trace(GN1)-3) - d*np.log(np.sqrt(I3N1)) + c/2*(np.sqrt(I3N1)-1)**2)*detJ*gaussWeight[k]
            
            # Derivative of the strain energy function
            Sigma_C = a*I                       # S_C = 2*DW/DC
            Sigma_G = b*I                       # S_G = 2*DW/DG
            CTest = np.linalg.norm(CN1-CN,'fro')
            if  CTest < 10**(-11):
                Sigma_I = -d/(2*I3N05)+c/2*(1-1/np.sqrt(I3N05))   # S_I =   2*DW/DI
                Sigma_I_I = d/(2*I3N05**2)+c/(2*(I3N05)**(3/2))    # S_I_I = 4*D2W/DI3^2
            else: #Greenspan formula '84
                Sigma_I = (-d*np.log(np.sqrt(I3N1)) + c/2*(np.sqrt(I3N1)-1)**2 - ( - d*np.log(np.sqrt(I3N)) + c/2*(np.sqrt(I3N)-1)**2))/(I3N1-I3N)
                Sigma_I_I = ((d/I3N1 - (c*(I3N1**(1/2) - 1))/I3N1**(1/2))/(I3N - I3N1) + (c*(I3N**(1/2) - 1)**2 - c*(I3N1**(1/2) - 1)**2 - 2*d*np.log(I3N**(1/2)) + 2*d*np.log(I3N1**(1/2)))/(I3N - I3N1)**2)
            
            # Second Piola Kirchhoff stress tensor
            SN05 = 2*(Sigma_C + wedgeAB(Sigma_G,CN05) + Sigma_I*1/3*(wedgeA(CN05) + GN05))
            SN05_v = np.array([SN05[0,0], SN05[1,1], SN05[2,2], SN05[0,1], SN05[1,2], SN05[0,2]])
            
            # Residual
            Re = Re + 1/2*BN05.T @ SN05_v*detJ*gaussWeight[k]
            
            # Tangent
            
            # Derivative of wedge(Sigma_G,CN05)
            D = (Sigma_G)
            SecDiffOperator = np.array([[0, D[2,2],  D[1,1],    0,        -D[2,1],  0],
                                        [D[2,2], 0,  D[0,0],    0,        0,       -D[2,0]],
                                        [D[1,1], D[0,0],  0, -D[1,0], 	  0,       0],
                                        [0,     0,   -D[1,0], -0.5*D[2,2], 0.5*D[2,0],  0.5*D[1,2]],
                                        [-D[2,1], 0,    0,    0.5*D[2,0],  -0.5*D[0,0],  0.5*D[0,1]],
                                        [0, -D[2,0],	0,    0.5*D[2,1],   0.5*D[1,0],  -0.5*D[1,1]]])
            # SecDiffOperator[4:6,4:6] = 0.5*SecDiffOperator(4:6,4:6);
            Kmat1 = SecDiffOperator
            
            # Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I
            D = (2/3*(CN05+0.5*CN1))
            SecDiffOperator = np.array([[0, D[2,2],  D[1,1],    0,        -D[2,1],  0],
                                        [D[2,2], 0,  D[0,0],    0,        0,       -D[2,0]],
                                        [D[1,1], D[0,0],  0, -D[1,0], 	  0,       0],
                                        [0,     0,   -D[1,0], -0.5*D[2,2], 0.5*D[2,0],  0.5*D[1,2]],
                                        [-D[2,1], 0,    0,    0.5*D[2,0],  -0.5*D[0,0],  0.5*D[0,1]],
                                        [0, -D[2,0],	0,    0.5*D[2,1],   0.5*D[1,0],  -0.5*D[1,1]]])
            # SecDiffOperator(4:6,4:6) = 0.5*SecDiffOperator(4:6,4:6);
            Kmat2 = Sigma_I*SecDiffOperator
            
            # Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
            GN05_v = np.array([GN05[0,0], GN05[1,1], GN05[2,2], GN05[0,1], GN05[1,2], GN05[0,2]])
            GMid_v = np.array([GMid[0,0], GMid[1,1], GMid[2,2], GMid[0,1], GMid[1,2], GMid[0,2]])
            GN1_v = np.array([GN1[0,0], GN1[1,1], GN1[2,2], GN1[0,1], GN1[1,2], GN1[0,2]])
            Kmat3 = Sigma_I_I*np.tensordot((GN05_v+2*GMid_v),GN1_v.T,axes=0)*1/3
            
            # Assembly of elasticity tensor
            ELA = Kmat1 + Kmat2 + Kmat3
            
            A1 = dNx.T@SN05@dNx*detJ*gaussWeight[k]
            MAT = np.zeros((numberOfDOFs,numberOfDOFs))
            for g in np.arange(0,dimension):
                MAT[np.ix_(range(g,numberOfDOFs,dimension),range(g,numberOfDOFs,dimension))] = A1

            Ke = Ke + 0.5*BN05.T@ELA@BN1*detJ*gaussWeight[k] + 0.5*MAT

        out['Re'][e] = Re
        out['pK'][e] = Ke.T.reshape(-1)
        out['ePot'][e] = ePot
        
    return out
