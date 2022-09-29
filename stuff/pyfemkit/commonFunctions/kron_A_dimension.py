import numpy as np
def kron_A_dimension(A, dimension):  # Simulates np.kron(A, np.eye(dimension))
    m,n = A.shape
    out = np.zeros((m,dimension,n,dimension),dtype=A.dtype)
    r = np.arange(dimension)
    out[:,r,:,r] = A
    out.shape = (m*dimension,n*dimension)
    return out