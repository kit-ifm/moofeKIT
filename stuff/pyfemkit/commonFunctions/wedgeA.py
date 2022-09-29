import numpy as np
def wedgeA(A):
# Creates tensor cross product.
     H = np.zeros((A.shape[0],A.shape[1]))
     H[0,:] = 2*np.cross(A[1,:],A[2,:])
     H[1,:] = 2*np.cross(A[2,:],A[0,:])
     H[2,:] = 2*np.cross(A[0,:],A[1,:])
     return H