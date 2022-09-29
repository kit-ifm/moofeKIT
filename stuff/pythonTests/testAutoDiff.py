import auto_diff
import numpy as np


# Residuum
def R(x,E,ddd):
    F_int = np.zeros((2,1))
    F_int[0,0] = x[0,0]**2 + x[1,0]**2 - E
    F_int[1,0] = x[0,0]**3 - x[1,0]**2 + 1
    
    F_ext = np.array([[1], [1]])
    
    out = (F_int - F_ext)*ddd
    return out

# Startwert
x1 = np.array([[1], [1]])

# Newton-Verfahren
newtonStep = 0
newtonMaximumSteps = 40
x = x1
E = 5
ddd = 27.0
while newtonStep <= newtonMaximumSteps:
    newtonStep = newtonStep + 1
    Res = R(x,E,ddd)
    normR = np.linalg.norm(Res)
    print('Iteration', newtonStep, ': norm(R)=', normR)
    if normR < 1e-8:
        break
    else:
        # Tangente
        with auto_diff.AutoDiff(x) as z:
            R_eval = R(z, E, ddd)
            K_T = auto_diff.jacobian(R_eval)
            
        # Inkrement
        dx = np.linalg.solve(K_T,Res)
        
        # Update
        x = x - dx

print('\nLösungsvektor x:\n', x)