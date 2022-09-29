import auto_diff
import numpy as np
from commonFunctions.expandEdof import expandEdof
from commonFunctions.computeDeterminant import computeDeterminant

def displacementSCMooneyRivlinAutoDiffDiscreteGradient(obj,setupObject):

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
        q = np.array([edN1.T.reshape(-1)]).T
        edN = qN[edof[e,:],].T
        # J = qR[np._ix(edof[e,:],range(dimension))].T @ dNrAll
        J = qR[edof[e,:],].T @ dNrAll
        # Run through all Gauss points
        for k in range(numberOfGausspoints):
            indx = range(dimension * k, dimension*(k+1))
            # detJ = np.linalg.det(J[:,indx].T)
            detJ = computeDeterminant(J[:,indx].T)
            if (detJ < 10*np.spacing(1)):
                print('Jacobi determinant equal or less than zero.'),  exit()
            checkIf = 'initial'
            re,checkIf = residualFunction(q,edN,J,dNrAll,indx,numberOfDOFs,dimension,a,b,c,d,I,checkIf)
            with auto_diff.AutoDiff(q) as z:
                re_eval,checkOut = residualFunction(z,edN,J,dNrAll,indx,numberOfDOFs,dimension,a,b,c,d,I,checkIf)
                ke = auto_diff.jacobian(re_eval)
            Re = Re + re.flatten()*detJ*gaussWeight[k]
            Ke = Ke + ke*detJ*gaussWeight[k]
        out['Re'][e] = Re
        out['pK'][e] = Ke.T.reshape(-1)
        out['ePot'][e] = ePot
    return out

def residualFunction(q,edN,J,dNrAll,indx,numberOfDOFs,dimension,a,b,c,d,I,checkIf):
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
    edN05 = 0.5*(edN+edN1)
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
    # ePot = ePot + (a*(np.trace(CN1)-3) + b*(np.trace(GN1)-3) - d*np.log(np.sqrt(I3N1)) + c/2*(np.sqrt(I3N1)-1)**2)*detJ*gaussWeight[k]
                
    # Derivative of the strain energy function
    Sigma_C = a*I                       # S_C = 2*DW/DC
    Sigma_G = b*I                       # S_G = 2*DW/DG
    # CTest = np.linalg.norm(CN1-CN,'fro')
    dC = CN1-CN
    # Sigma_I = -d/(2*I3N05)+c/2*(1-1/np.sqrt(I3N05))   # S_I =   2*DW/DI
    CTest = np.sqrt(sum(sum(dC**2)))
    if (checkIf == 'initial'):
        if  CTest < 10**(-11):
            checkIf = True
        else: #Greenspan formula '84
            checkIf = False
    if checkIf:
        Sigma_I = -d/(2*I3N05)+c/2*(1-1/np.sqrt(I3N05))   # S_I =   2*DW/DI
    else:
        Sigma_I = (-d*np.log(np.sqrt(I3N1)) + c/2*(np.sqrt(I3N1)-1)**2 - ( - d*np.log(np.sqrt(I3N)) + c/2*(np.sqrt(I3N)-1)**2))/(I3N1-I3N)

    # if  CTest < 10**(-11):
    #      Sigma_I = -d/(2*I3N05)+c/2*(1-1/np.sqrt(I3N05))   # S_I =   2*DW/DI
    # else: #Greenspan formula '84
    #     Sigma_I = (-d*np.log(np.sqrt(I3N1)) + c/2*(np.sqrt(I3N1)-1)**2 - ( - d*np.log(np.sqrt(I3N)) + c/2*(np.sqrt(I3N)-1)**2))/(I3N1-I3N)

                
    # Second Piola Kirchhoff stress tensor
    SN05 = 2*(Sigma_C + wedgeAB(Sigma_G,CN05) + Sigma_I*1/3*(wedgeA(CN05) + GN05))
    SN05_v = np.array([[SN05[0,0]], [SN05[1,1]], [SN05[2,2]], [SN05[0,1]], [SN05[1,2]], [SN05[0,2]]])
                
    # Residual
    re = 0.5*BN05.T@SN05_v
    
    return re,checkIf

def wedgeA(A):
# %% Creates tensor cross product.
#     %
#     % Syntax
#     % A,B vector or matrix
#     %
#     % [H] = wedge(FAkt,FAkt)      
#     %
#     % Description
#     %
#     % Tensor cross product.
#     % 16.01.2015 C. Bilgen

     H = np.zeros((A.shape[0],A.shape[1]))
    
     H[0,0] = 2*(A[1,1]*A[2,2]-A[1,2]*A[2,1])
     H[0,1] = 2*(A[1,2]*A[2,0]-A[1,0]*A[2,2])
     H[0,2] = 2*(A[1,0]*A[2,1]-A[1,1]*A[2,0])
     
     H[1,0] = 2*(A[2,1]*A[0,2]-A[2,2]*A[0,1])
     H[1,1] = 2*(A[2,2]*A[0,0]-A[2,0]*A[0,2])
     H[1,2] = 2*(A[2,0]*A[0,1]-A[2,1]*A[0,0])

     H[2,0] = 2*(A[0,1]*A[1,2]-A[0,2]*A[1,1])
     H[2,1] = 2*(A[0,2]*A[1,0]-A[0,0]*A[1,2])
     H[2,2] = 2*(A[0,0]*A[1,1]-A[0,1]*A[1,0])

     return H
    
    # TODO
    # if size(A) == size(B) 
    #     if  A == B
    #         H = zeros(size(A));
    #         H(1,:) = 2*cross(A(2,:),A(3,:));
    #         H(2,:) = 2*cross(A(3,:),A(1,:));
    #         H(3,:) = 2*cross(A(1,:),A(2,:));
    #     else
    #         H = zeros(size(A));
    #         H(1,:) = cross(A(2,:),B(3,:)) + cross(B(2,:),A(3,:));
    #         H(2,:) = cross(B(3,:),A(1,:)) + cross(A(3,:),B(1,:));
    #         H(3,:) = cross(A(1,:),B(2,:)) + cross(B(1,:),A(2,:));
    #     end
    # elseif isvector(A) || isvector(B)
    #     if isvector(A) == 1 && isvector(B) == 0
    #         H = zeros(size(B));
    #         H(:,1) = cross(A,B(:,1));
    #         H(:,2) = cross(A,B(:,2));
    #         H(:,3) = cross(A,B(:,3));
    #     elseif isvector(A) == 0 && isvector(B) == 1
    #         H = zeros(size(A));
    #         H(1,:) = cross(A(1,:),B);       
    #         H(2,:) = cross(A(2,:),B);
    #         H(3,:) = cross(A(3,:),B);
    #     else 
    #         error('matrix dimensions do not agree!')
    #     end
       
    # else
    #     error('Matrix dimensions do not agree!')
    # end
# end

def wedgeAB(A,B):
# %% Creates tensor cross product.
#     %
#     % Syntax
#     % A,B vector or matrix
#     %
#     % [H] = wedge(FAkt,FAkt)      
#     %
#     % Description
#     %
#     % Tensor cross product.
#     % 16.01.2015 C. Bilgen

     H = np.zeros((A.shape[0],A.shape[1]))
    
     H[0,0] = (A[1,1]*B[2,2]-A[1,2]*B[2,1]) + (B[1,1]*A[2,2]-B[1,2]*A[2,1])
     H[0,1] = (A[1,2]*B[2,0]-A[1,0]*B[2,2]) + (B[1,2]*A[2,0]-B[1,0]*A[2,2])
     H[0,2] = (A[1,0]*B[2,1]-A[1,1]*B[2,0]) + (B[1,0]*A[2,1]-B[1,1]*A[2,0])
        
     H[1,0] = (B[2,1]*A[0,2]-B[2,2]*A[0,1]) + (A[2,1]*B[0,2]-A[2,2]*B[0,1])
     H[1,1] = (B[2,2]*A[0,0]-B[2,0]*A[0,2]) + (A[2,2]*B[0,0]-A[2,0]*B[0,2])
     H[1,2] = (B[2,0]*A[0,1]-B[2,1]*A[0,0]) + (A[2,0]*B[0,1]-A[2,1]*B[0,0])

     H[2,0] = (A[0,1]*B[1,2]-A[0,2]*B[1,1]) + (B[0,1]*A[1,2]-B[0,2]*A[1,1])
     H[2,1] = (A[0,2]*B[1,0]-A[0,0]*B[1,2]) + (B[0,2]*A[1,0]-B[0,0]*A[1,2])
     H[2,2] = (A[0,0]*B[1,1]-A[0,1]*B[1,0]) + (B[0,0]*A[1,1]-B[0,1]*A[1,0])

     return H
