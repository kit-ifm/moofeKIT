import numpy as np
def kron_dimension_A(A, dimension):  # Simulates np.kron(np.eye(dimension), A)
    m,n = A.shape
    out = np.zeros((dimension,m,dimension,n),dtype=A.dtype)
    r = np.arange(dimension)
    out[r,:,r,:] = A
    out.shape = (m*dimension,n*dimension)
    return out
