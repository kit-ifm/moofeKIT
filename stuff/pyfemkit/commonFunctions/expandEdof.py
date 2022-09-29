import numpy as np
def expandEdof(edof):
    numberOfElements = len(edof)
    # edof1 = (np.kron(np.ones(numberOfElements),edof)).astype('i')
    # edof2 = (np.kron(edof,np.ones(numberOfElements))).astype('i')
    edof1 = np.tile(edof,numberOfElements)
    edof2 = np.repeat(edof,numberOfElements)
    return edof1, edof2