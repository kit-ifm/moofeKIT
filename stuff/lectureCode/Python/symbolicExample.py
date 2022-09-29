# Tutorial: https://docs.sympy.org/latest/tutorial/index.html

from sympy import *


# SUBFUNCTIONS
def GradientOfVec(u):
    X1, X2, X3 = symbols('X1 X2 X3')
    X = Matrix([[X1, X2, X3]])
    out = zeros(3, 3)
    for i in range(3):
        out = out + zeros(3, i).row_join(diff(u, X[i])).row_join(zeros(3, 2-i))
    return out


def DivergenceOfTensor(A):
    X1, X2, X3 = symbols('X1 X2 X3')
    X = Matrix([[X1, X2, X3]])
    out = zeros(3, 1)
    I = eye(3)
    for i in range(3):
        out = out + diff(A, X[i]) * I[:, i]
    return out


# MAIN
X1, X2, X3 = symbols('X1 X2 X3')

# displacement (analytic)
xSymbolic = Matrix([[X1+10*X1**2], [X2], [X3]])

Fana = GradientOfVec(xSymbolic)

# constitutive equations
# Psi = a*(Trace(F.T*F)-3)
DPsi_F = 12*2*Fana
Pana = DPsi_F

# body force
BSymbolic = -DivergenceOfTensor(Pana)

# substitution
FanaS = Fana.subs([(X1, 10), (X2, 2), (X3, 0)])
