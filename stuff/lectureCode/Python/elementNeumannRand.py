import numpy as np

def elementNeumannRand(ed_rand, ir, RB_2):
    # Lastroutine (Neumannrand)
    #
    #   siehe Elementroutine


    # Anzahl der Integrationspunkte
    ngp = ir

    # Fallunterscheidung - Stuetzpunkte, Gewichte
    if ir == 1:
        g1 = 0.0
        w1 = 2.0
        gp = g1
        w = w1

    elif ir == 2:
        g1 = 0.577350269189626
        w1 = 1
        gp = np.array([-g1, g1])
        w  = np.array([w1, w1])

    elif ir == 3:
        g1 = 0.774596669241483
        g2 = 0
        w1 = 0.555555555555555
        w2 = 0.888888888888888
        gp = np.array([-g1, g2, g1])
        w  = np.array([w1, w2, w1])

    else:
        print('Used number of intergratins poits not implemented')
        return

    xsi = gp[:]

    ## Lineare Formfunktionen
    #   N = [   N1(xi1)     N2(xi1)
    #           N1(xi2)     N2(xi2) ]
    N = np.zeros((2,2))
    N[:, 0] = (1 - xsi) / 2
    N[:, 1] = (1 + xsi) / 2

    # Elementlaenge
    Le = np.linalg.norm(np.array([ed_rand[2], ed_rand[3]]) - np.array([ed_rand[0], ed_rand[1]]))

    # Jacobimatrix/-determinante
    J = Le / 2
    detJ = J

    ## Elementlastvektor aufbauen
    Fe = np.zeros(4)
    for ii in np.arange(0,ngp):
        # Schleife ueber alle Gausspunkte   
        Fe[np.array([0,2])] = Fe[np.array([0,2])] + N[ii, :].T * RB_2[0] * detJ * w[ii]
        Fe[np.array([1,3])] = Fe[np.array([1,3])] + N[ii, :].T * RB_2[1] * detJ * w[ii]
        
    return np.atleast_2d(Fe)