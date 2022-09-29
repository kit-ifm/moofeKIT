import numpy as np
#from scipy.sparse.linalg import spsolve
from scikits.umfpack import spsolve, splu

def solveq(K, F, bc):
    ## Loesen des lineare Gleichungssystems unter Beruecksichtiung Dirichlet-Randbedingungen
    #
    # q = solveq(K,F,bc)
    #
    # q: Loesungsvektor
    # K: Steifigkeitsmatrix
    # F: Lastvektor
    # bc: Dirichlet-Freiheitsgrade
    #
    # bc = [idxDiri1 valueDiri1
    #       idxDiri2 valueDiri2
    #            :
    #       idxDiriN valueDiriN]

    if bc.size != 0:
        # Eingaben einlesen
        dofDiri = (bc[:, 0] - np.ones(bc[:, 0].shape)).astype('i')
        qDiri = bc[:, 1]

        # Freiheitsgrade mit Dirichlet-Raender identifizieren
        nDof = np.size(F)
        DOSOLVE = np.full(nDof, True)
        DOSOLVE[dofDiri] = False

        # Nach Freiheitsgraden, die nicht Dirichlet-Rand sind loesen, bei
        # anderen Werten Dirichletrand einsetzen
        q = np.zeros(np.shape(K)[0])
        q[DOSOLVE] = spsolve(K[np.ix_(DOSOLVE, DOSOLVE)], F[DOSOLVE] - np.atleast_2d(K[np.ix_(DOSOLVE, np.logical_not(DOSOLVE))] @ qDiri).T)
        q[np.logical_not(DOSOLVE)] = qDiri
    else:
        # Nach allen Werten loesen
        q = np.linalg.solve(K, F)

    return q
