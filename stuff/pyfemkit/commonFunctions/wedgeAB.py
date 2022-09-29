import numpy as np
def wedgeAB(A,B):
# Creates tensor cross product.
     H = np.zeros((A.shape[0],A.shape[1]))
     H[0,:] = np.cross(A[1,:],B[2,:]) + np.cross(B[1,:],A[2,:])
     H[1,:] = np.cross(B[2,:],A[0,:]) + np.cross(A[2,:],B[0,:])
     H[2,:] = np.cross(A[0,:],B[1,:]) + np.cross(B[0,:],A[1,:])
     return H
