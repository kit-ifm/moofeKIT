import numpy as np

def elementSRFull(ed,ir,C):
    # Elementroutine
    #
    #   Eingang: ed - Zeile der Elementgeometriedatentabelle
    #               ed = [ x1 y1 x2 y2 x3 y3 x4 y4 ];
    #            ir - Parameter Integrationsstuetzpunkte
    #            C  - Materialmatrix


    # Anzahl der Integrationspunkte
    ngp = ir*ir

    # Fallunterscheidung - Stuetzpunkte, Gewichte
    if ir == 1:
        g1 = 0.0 
        w1 = 2.0
        gp = np.array([g1, g1])
        w  = np.array([w1, w1])
    
    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([[-g1, -g1], [g1, -g1], [-g1, g1], [g1, g1]])
        w  = np.array([[w1, w1], [w1, w1], [w1, w1], [w1, w1]])
    
    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0
        w1 = 0.555555555555555 
        w2 = 0.888888888888888
        gp = np.array([[-g1, -g1], [-g2, -g1], [g1, -g1], [-g1, g2], [g2, g2], [g1, g2], [-g1, g1], [g2, g1], [g1, g1]])
        w  = np.array([[w1, w1], [w2, w1], [w1, w1], [w1, w2], [w2, w2], [w1, w2], [w1, w1], [w2, w1], [w1, w1]])
    
    else:
        print('Used number of integration points not implemented');
        return
    
    wp  = w[:,0]*w[:,1]
    
    xsi  = gp[:,0]
    eta = gp[:,1]
    r2  = ngp*2
    
    ## Formfunktionen 
    # N = [ N1(xi1,eta1) N2(xi1,eta1) N3(xi1,eta1) N4(xi1,eta1) 
    #       N1(xi2,eta2) N2(xi2,eta2) N3(xi2,eta2) N4(xi2,eta2) 
    #       N1(xi3,eta3) N2(xi3,eta3) N3(xi3,eta3) N4(xi3,eta3) 
    #       N1(xi4,eta4) N2(xi4,eta4) N3(xi4,eta4) N4(xi4,eta4) ]
    #
    #   Tauchen hier nicht in der STEMA auf, deswegen auskommentiert.
    
    # N = np.zeros((4,4))
    # % N(:,1) = (1-xi).*(1-eta)/4;  
    # % N(:,2) = (1+xi).*(1-eta)/4;
    # % N(:,3) = (1+xi).*(1+eta)/4;  
    # % N(:,4) = (1-xi).*(1+eta)/4;     
    
    
    ## Ableitungen der Formfunktionen
    # dNr = [ dN1dxi (eta1)   dN2dxi (eta1)    dN3dxi (eta1)     dN4dxi(eta1)
    #         dN1deta(xi1)    dN2deta(xi1)     dN3deta(xi1)      dN4deta(xi1)
    #         dN1dxi (eta2)   dN2dxi (eta2)    dN3dxi (eta2)     dN4dxi(eta2)
    #         dN1deta(xi2)    dN2deta(xi2)     dN3deta(xi2)      dN4deta(xi2)
    #         dN1dxi (eta3)   dN2dxi (eta3)    dN3dxi (eta3)     dN4dxi(eta3)
    #         dN1deta(xi3)    dN2deta(xi3)     dN3deta(xi3)      dN4deta(xi3)
    #         dN1dxi (eta4)   dN2dxi (eta4)    dN3dxi (eta4)     dN4dxi(eta4)
    #         dN1deta(xi4)    dN2deta(xi4)     dN3deta(xi4)      dN4deta(xi4) ]
    dNr = np.zeros((8,4))

    # Partielle Ableitungen der Formfunktionen nach xi
    dNr[0:r2:2,0] = -(1-eta)/4;
    dNr[0:r2:2,1] =  (1-eta)/4;
    dNr[0:r2:2,2] =  (1+eta)/4;     
    dNr[0:r2:2,3] = -(1+eta)/4;

    # Partielle Ableitungen der Formfunktionen nach eta
    dNr[1:r2+1:2,0] = -(1-xsi)/4;   
    dNr[1:r2+1:2,1] = -(1+xsi)/4;
    dNr[1:r2+1:2,2] =  (1+xsi)/4;   
    dNr[1:r2+1:2,3] =  (1-xsi)/4;


    ## Jacobimatrizen (fuer alle Gausspunkte)
    #   J = [ J(xi1,eta1) J(xi2,eta2) J(xi3,eta3) J(xi4,eta4) ]
    #J = np.array([ed[0:2:ed.shape[0]-1], ed[1:2:ed.shape[0]]]) @ dNr.T
    J = np.array([ed[0:ed.shape[0]-5:2], ed[1:ed.shape[0]-4:2]]) @ dNr.T
    
    eta0=0
    xsi0=0
    dNr0 = np.zeros((8,4))
    
     # Partielle Ableitungen der Formfunktionen nach xi
    dNr0[0:r2:2,0] = -(1-eta0)/4;
    dNr0[0:r2:2,1] =  (1-eta0)/4;
    dNr0[0:r2:2,2] =  (1+eta0)/4;     
    dNr0[0:r2:2,3] = -(1+eta0)/4;

    # Partielle Ableitungen der Formfunktionen nach eta
    dNr0[1:r2+1:2,0] = -(1-xsi0)/4;   
    dNr0[1:r2+1:2,1] = -(1+xsi0)/4;
    dNr0[1:r2+1:2,2] =  (1+xsi0)/4;   
    dNr0[1:r2+1:2,3] =  (1-xsi0)/4;
     ## Jacobimatrizen (fuer alle Gausspunkte)
    #   J = [ J(xi1,eta1) J(xi2,eta2) J(xi3,eta3) J(xi4,eta4) ]
    #J0h = np.array([ed[0:ed.shape[0]-1:2], ed[1:ed.shape[0]:2]])
    J0 = np.array([ed[0:ed.shape[0]-5:2], ed[1:ed.shape[0]-4:2]]) @ dNr0.T
    #J0 = np.array([ed[0:2:ed.shape[0]-1], ed[1:2:ed.shape[0]]]) @ dNr.T
    
    ## Elementsteifigkeitmatrizen/-lastvektoren aufbauen
    # Allozieren
    Ke = np.zeros((8,8))
    He = np.zeros((4,4))
    Te = np.zeros((4,8))
    Fe = np.zeros(12)
    
    J011 = J0[0,0]
    J012 = J0[0,1]
    J021 = J0[1,0]
    J022 = J0[1,1]
    F0 = np.array([[J011*J011,J021*J012,2*J011*J012],[J012*J021,J022*J022,2*J021*J022],[J011*J021,J012*J022,J011*J022+J012*J021]])
    
    j0=np.linalg.det(J0[:,[0,1]])
    
    for i in np.arange(0,ngp):   
        # Schleife ueber alle Gausspunkte
        
        # Zugriffsindex fuer korrekte Formfunktionen/Jacobimatrix
        indx    = np.array([ 2*i, 2*i+1 ])
        
        # Jacobimatrix/-determinante
        JT      = J[:,indx].T
        detJ    = np.linalg.det(JT)
    
        # Plausibilitaetstest
        if detJ < 10*np.finfo(float).eps:
            print('Jacobideterminant equal or less than zero!')
    
        # Ableitung Formfunktion nach globalen Koordinaten
        #   dNdx = [    dN1dx   dN2dx   dN3dx   dN4dx
        #               dN1dy   dN2dy   dN3dy   dN4dy   ]
        dNx = np.linalg.solve(JT,dNr[indx,:])
        
    
        # Knotenoperatormatrix
        # B = [   dN1dx   0       dN2dx   0       dN3dx   0       dN4dx   0
        #         0       dN1dy   0       dN2dy   0       dN3dy   0       dN4dy
        #         dN1dy   dN1dx   dN2dy   dN2dx   dN3dy   dN3dx   dN4dy   dN4dx ]    
        B = np.zeros((3,8))
        B[0,0:7:2] = dNx[0,:]
        B[1,1:8:2] = dNx[1,:]
        B[2,0:7:2] = dNx[1,:]
        B[2,1:8:2] = dNx[0,:]
        
        
        # Elementsteifigkeitsmatrix
        Ke = Ke + B.T @ C @ B * detJ * wp[i]
        
        E=np.array([[xsi[i],0,0,0],[0,eta[i],0,0],[0,0,xsi[i],eta[i]]])
        G=j0/detJ * np.linalg.inv(F0.T) @ E
        He = He + G.T @ C @ G * detJ * wp[i]
        Te = Te + G.T @ C @ B * detJ * wp[i]
    KeFull1 = np.concatenate((Ke,(Te.T)),axis=1)
    KeFull2 = np.concatenate((Te,He),axis=1)
    KeFull=np.concatenate((KeFull1,KeFull2))
        
    return KeFull, Fe