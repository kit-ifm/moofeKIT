import numpy as np

# from setupClass import setupClass
# from dofClass import dofClass

def lagrangeShapeFunctions(numNodes, numberOfGausspoints, dimension):

    out = {"N": None,
           "dNr": None,
           "dNr2": None,
           "gaussWeight": None,
           "gaussPoint": None
           }
    g1 = 0.577350269189626
    g2 = 0.774596669241483
    g3 = 0.339981043584856
    g4 = 0.861136311594053
         
    w1 = 0.555555555555556
    w2 = 0.888888888888889
    w3 = 0.652145154862546
    w4 = 0.347854845137454

    if (dimension==1):
        #Gausspoints
        if (numberOfGausspoints == 1):
            gaussPoints = 0
            gaussWeight = 2
        elif (numberOfGausspoints == 2):
            g = 1/np.sqrt(3)
            gaussPoints =  np.array([-g,g])
            gaussWeight =  np.array([1,1])
        elif (numberOfGaussPoints == 3):
            g = np.sqrt(3/5)
            gaussPoints =  np.array([-g,0,g])
            gaussWeight =  np.array([5/9,8/9,5/9])
        else:
            print('number of Gausspoints not implemented')

        #Ansatzfunctions at gausspoints
        xsi = gaussPoint[0,:]
        N = np.zeros((numberOfGausspoints,numNodes))
        dNr = N
        dNr2 = N

        if (numNodes==1):
            N[:,0] = 1
            dNr[:,0] = 0
            dNr[:,0] = 0
        elif (numNodes==2):
            #shape functions
            N[:,0] = 1/2 * (1-xsi)
            N[:,1] = 1/2 * (1+xsi)
            #derivative wrt xsi
            dNr[:,0] = xsi - 1/2
            dNr[:,1] = xsi + 1/2
            dNr[:,2] = -2 *xsi
            #second derivative wrt xsi
            dNr2[:,0] = 1
            dNr2[:,1] = 1
            dNr2[:,2] = -2
        else:
            print('Order of ansatzfunctions not implemented')
    
    elif (dimension == 2):
        # Gausspoints
        if (numberOfGausspoints == 1):
            gaussPoints =  np.array([[0],[0]])
            gaussWeight = 2
        elif (numberOfGausspoints == 4):
            g = 1/np.sqrt(3)
            gaussPoints =  np.array([[-g,g,g,-g],[-g,-g,g,g]])
            gaussWeight =  np.array([1,1,1,1])
        elif (numberOfGausspoints == 9):
            g = np.sqrt(3/5)
            gaussPoints =  np.array([[-g,0,g,-g,0,g,-g,0,g],[-g,-g,-g,0,0,0,g,g,g]])
            gaussWeight =  np.array([25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81])
        else:
            print('number of Gausspoints not implemented')

        #Ansatzfunctions at gausspoints
            xsi = gaussPoints[0,:]
            eta = gaussPoints[1,:]
            r2 = numberOfGausspoints * 2
            r3 = numberOfGausspoints * 3
            N = np.zeros((numberOfGausspoints,numNodes))
            dNr = np.zeros((r2,numNodes))
            dNr2 = np.zeros((r3,numNodes))

            if (numNodes == 1):
                N[:,:] = 1
                dNr[:,:] = 0
                dNr2[:,:] = 0
            elif (numNodes == 4):
                #shape functions
                N[:,0] = (1-xsi)*(1-eta)/4
                N[:,1] = (1+xsi)*(1-eta)/4
                N[:,2] = (1+xsi)*(1+eta)/4
                N[:,3] = (1-xsi)*(1+eta)/4
                #derivative wrt xsi
                dNr[0:r2:2,0] = -(1-eta)/4
                dNr[0:r2:2,1] = (1-eta)/4
                dNr[0:r2:2,2] = (1+eta)/4
                dNr[0:r2:2,3] = (1+eta)/4
                #derivative wrt eta
                dNr2[1:r2:2,0] = -(1-xsi)/4
                dNr2[1:r2:2,1] = -(1+xsi)/4
                dNr2[1:r2:2,2] = (1+xsi)/4
                dNr2[1:r2:2,3] = (1-xsi)/4
            elif (numNodes == 9):
                #shape functions
                N[:,0] = 1/4 * (xsi**2 - xsi) * (eta**2 - eta)
                N[:,1] = 1/4 * (xsi**2 + xsi) * (eta**2 - eta)
                N[:,2] = 1/4 * (xsi**2 + xsi) * (eta**2 + eta)
                N[:,3] = 1/4 * (xsi**2 - xsi) * (eta**2 + eta)
                N[:,4] = 1/2 * (1 - xsi**2) * (eta**2),(eta-1)
                N[:,5] = 1/2 * xsi * (xsi+1) * (1-eta**2)
                N[:,6] = 1/2 * (1-xsi**2) * eta * (eta+1)
                N[:,7] = 1/2 * xsi * (xsi-1) * (1-eta**2)
                N[:,8] = 1/2 * (1-xsi**2) * (1-eta**2)
                #derivative wrt xsi
                dNr[0:r2:2,0] = 1/4 * (2*xsi - 1) * (eta**2 - eta)
                dNr[0:r2:2,1] = 1/4 * (2*xsi + 1) * (eta**2 - eta)
                dNr[0:r2:2,2] = 1/4 * (2*xsi + 1) * (eta**2 + eta)
                dNr[0:r2:2,3] = 1/4 * (2*xsi - 1) * (eta**2 + eta)
                dNr[0:r2:2,4] = -xsi * eta * (eta-1) 
                dNr[0:r2:2,5] = 1/2 * (2*xsi + 1) * (1 - eta**2)
                dNr[0:r2:2,6] = -xsi * eta * (eta + 1)
                dNr[0:r2:2,7] = 1/2 * (2*xsi - 1) * (1 - eta**2)
                dNr[0:r2:2,8] = -2 * xsi * (1 - eta**2)
                #derivative wrt eta
                dNr[1:r2:2,0] = 1/4 * (xsi**2 - xsi) * (2*eta - 1)
                dNr[1:r2:2,1] = 1/4 * (xsi**2 + xsi) * (2*eta - 1)
                dNr[1:r2:2,2] = 1/4 * (xsi**2 + xsi) * (2*eta + 1)
                dNr[1:r2:2,3] = 1/4 * (xsi**2 - xsi) * (2*eta + 1)
                dNr[1:r2:2,4] = 1/2 * (1-xsi**2) * (2*eta - 1)
                dNr[1:r2:2,5] = -xsi * (xsi+1) * eta
                dNr[1:r2:2,6] = 1/2 * (1-xsi**2) * (2*eta + 1)
                dNr[1:r2:2,7] = -xsi * (xsi - 1) * eta
                dNr[1:r2:2,8] = -2 * (1-xsi**2) * eta
            else:
                print('Order of ansatzfunctions not implemented')
    elif (dimension == 3):
        if (numberOfGausspoints == 1):
            gaussPoints = np.array([0,0,0])
            gaussWeight = 8
        elif (numberOfGausspoints == 8):
            gaussPoints =  np.array([[-g1,-g1,-g1],
                                     [g1,-g1,-g1],
                                     [-g1,g1,-g1],
                                     [g1,g1,-g1],
                                     [-g1,-g1,g1],
                                     [g1,-g1,g1],
                                     [-g1,g1,g1],
                                     [g1,g1,g1]])
            gaussWeight = np.ones((8,1))
        elif (numberOfGausspoints == 28):
            gaussPoints =  np.array([[-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2,-g2,0,g2],
                           [-g2,-g2,-g2,0,0,0,g2,g2,g2,-g2,-g2,-g2,0,0,0,g2,g2,g2,-g2,-g2,-g2,0,0,0,g2,g2,g2],
                           [-g2,-g2,-g2,-g2,-g2,-g2,-g2,-g2,-g2,0,0,0,0,0,0,0,0,0,g2,g2,g2,g2,g2,g2,g2,g2,g2]])
            gaussWeight =  np.array([[25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81].dT*w1,
                           [25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81].dT*w2,
                           [25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81].dT*w1])
        else:
            print('Number of Gauss points not implemented')

        #Support
        r3 = numberOfGausspoints * 3
        xsi = gaussPoints[:,0]
        eta = gaussPoints[:,1]
        zeta = gaussPoints[:,2]
        N = np.zeros((numberOfGausspoints,numNodes))
        dNr = np.zeros((r3,numNodes))
        dNr2 = np.zeros((r3,numNodes))

        if (numNodes == 8):
            #shape functions
            N[:,0] = (1-xsi)*(1-eta)*(1-zeta)/8
            N[:,1] = (1+xsi)*(1-eta)*(1-zeta)/8
            N[:,2] = (1+xsi)*(1+eta)*(1-zeta)/8
            N[:,3] = (1-xsi)*(1+eta)*(1-zeta)/8
            N[:,4] = (1-xsi)*(1-eta)*(1+zeta)/8
            N[:,5] = (1+xsi)*(1-eta)*(1+zeta)/8
            N[:,6] = (1+xsi)*(1+eta)*(1+zeta)/8
            N[:,7] = (1-xsi)*(1+eta)*(1+zeta)/8

            #derivative wrt xsi
            dNr[0:r3:3,0] = (-(1-eta)) * (1-zeta)/8
            dNr[0:r3:3,1] = (1-eta) * (1-zeta)/8
            dNr[0:r3:3,2] = (1+eta) * (1-zeta)/8
            dNr[0:r3:3,3] = (-(1+eta)) * (1-zeta)/8
            dNr[0:r3:3,4] = (-(1-eta)) * (1+zeta)/8
            dNr[0:r3:3,5] = (1-eta) * (1+zeta)/8
            dNr[0:r3:3,6] = (1+eta) * (1+zeta)/8
            dNr[0:r3:3,7] = (-(1+eta)) * (1+zeta)/8

            #derivative wrt eta
            dNr[1:r3:3,0] = (-(1-xsi)) * (1-zeta)/8
            dNr[1:r3:3,1] = (-(1+xsi)) * (1-zeta)/8
            dNr[1:r3:3,2] = (1+xsi) * (1-zeta)/8
            dNr[1:r3:3,3] = (1-xsi) * (1-zeta)/8
            dNr[1:r3:3,4] = (-(1-xsi)) * (1+zeta)/8
            dNr[1:r3:3,5] = (-(1+xsi)) * (1+zeta)/8
            dNr[1:r3:3,6] = (1+xsi) * (1+zeta)/8
            dNr[1:r3:3,7] = (1-xsi) * (1+zeta)/8

            #derivative wrt zeta
            dNr[2:r3:3,0] = (-(1-xsi)) * (1-eta)/8
            dNr[2:r3:3,1] = (-(1+xsi)) * (1-eta)/8
            dNr[2:r3:3,2] = (-(1+xsi)) * (1+eta)/8
            dNr[2:r3:3,3] = (-(1-xsi)) * (1+eta)/8
            dNr[2:r3:3,4] = (1-xsi) * (1-eta)/8
            dNr[2:r3:3,5] = (1+xsi) * (1-eta)/8
            dNr[2:r3:3,6] = (1+xsi) * (1+eta)/8
            dNr[2:r3:3,7] = (1-xsi) * (1+eta)/8

            dNr2[0:r3:3,:] = 0

            dNr2[1:r3:3,0] = (1-zeta)/8
            dNr2[1:r3:3,1] = -(1-zeta)/8
            dNr2[1:r3:3,2] = (1-zeta)/8
            dNr2[1:r3:3,3] = -(1-zeta)/8
            dNr2[1:r3:3,4] = (1+zeta)/8
            dNr2[1:r3:3,5] = -(1+zeta)/8
            dNr2[1:r3:3,6] = (1+zeta)/8
            dNr2[1:r3:3,7] = -(1+zeta)/8

            dNr2[2:r3:3,0] = (1-eta)/8
            dNr2[2:r3:3,1] = -(1-eta)/8
            dNr2[2:r3:3,2] = -(1+eta)/8
            dNr2[2:r3:3,3] = (1+eta)/8
            dNr2[2:r3:3,4] = -(1-eta)/8
            dNr2[2:r3:3,5] = (1-eta)/8
            dNr2[2:r3:3,6] = (1+eta)/8
            dNr2[2:r3:3,7] = -(1+eta)/8

            dNr2[1:r3:3,0] = (1-zeta)/8
            dNr2[1:r3:3,1] = -(1-zeta)/8
            dNr2[1:r3:3,2] = (1-zeta)/8
            dNr2[1:r3:3,3] = -(1-zeta)/8
            dNr2[1:r3:3,4] = (1+zeta)/8
            dNr2[1:r3:3,5] = -(1+zeta)/8
            dNr2[1:r3:3,6] = (1+zeta)/8
            dNr2[1:r3:3,7] = -(1+zeta)/8

            dNr2[1:r3:3,:] = 0

            dNr2[2:r3:3,0] = (1-xsi)/8
            dNr2[2:r3:3,1] = (1+xsi)/8
            dNr2[2:r3:3,2] = -(1+xsi)/8
            dNr2[2:r3:3,3] = -(1-xsi)/8
            dNr2[2:r3:3,4] = -(1-xsi)/8
            dNr2[2:r3:3,5] = -(1+xsi)/8
            dNr2[2:r3:3,6] = (1+xsi)/8
            dNr2[2:r3:3,7] = (1-xsi)/8

            dNr2[0:r3:3,0] = (1-eta)/8
            dNr2[0:r3:3,1] = -(1-eta)/8
            dNr2[0:r3:3,2] = -(1+eta)/8
            dNr2[0:r3:3,3] = (1+eta)/8
            dNr2[0:r3:3,4] = -(1-eta)/8
            dNr2[0:r3:3,5] = (1-eta)/8
            dNr2[0:r3:3,6] = (1+eta)/8
            dNr2[0:r3:3,7] = -(1+eta)/8
            dNr2[1:r3:3,0] = (1-xsi)/8
            dNr2[1:r3:3,1] = (1+xsi)/8
            dNr2[1:r3:3,2] = -(1+xsi)/8
            dNr2[1:r3:3,3] = -(1-xsi)/8
            dNr2[1:r3:3,4] = -(1-xsi)/8
            dNr2[1:r3:3,5] = -(1+xsi)/8
            dNr2[1:r3:3,6] = (1+xsi)/8
            dNr2[1:r3:3,7] = (1-xsi)/8

            dNr2[2:r3:3,:] = 0

        else:
            print('not implemented')

    else:
        print('dimension not implemented')

    #output
    out["gaussWeight"] = gaussWeight
    out["gaussPoint"] = gaussPoints
    out["N"] = N
    out["dNr"] = dNr
    out["dNr2"] = dNr2

    return out
