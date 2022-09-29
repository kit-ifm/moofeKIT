import numpy as np

def assem(O, Oe, dof):
    ## Finite element assembly
    # 
    # K = assem(K,Ke,dof)
    # F = assem(F,Fe,dof)
    #
    # K: Steifigkeitsmatrix
    # Ke: Elementsteifigkeitsmatrix
    #
    # F: Lastvektor
    # Fe: Elementlastvektor
    #
    # dof: Indices der Elementfreiheitsgrade
    #
    # Zusammenbau der Finite-Elemente-Matrizes oder FE-Vektoren

    dof = dof - np.ones(dof.shape,dtype="i")

    if Oe.ndim == 2 and Oe.shape[0] > 1:
        # Matrix-Assembly
        # K(dof,dof) = K(dof,dof) + Ke
        O[np.ix_(dof, dof)] = O[np.ix_(dof, dof)] + Oe
    elif Oe.ndim == 2:
        # Vektor-Assembly
        # F(dof) = F(dof) + Fe;
        O[dof,0] = O[dof,0].T + Oe
    else:
        # Vektor-Assembly
        # F(dof) = F(dof) + Fe;
        O[dof] = O[dof] + Oe
        
    return O