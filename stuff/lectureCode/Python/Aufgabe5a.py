## Importieren benoetigter Bibliotheken und Funktionen
# numPy & SciPy
import numpy as np
from scipy.sparse import csr_matrix

# Funktionen
from netz import netz
from extract import extract
from draw import draw
from element import element
from elementSRFull import elementSRFull
from assem import assem
from elementNeumannRand import elementNeumannRand
from solveq import solveq
from nachlauf import nachlauf
from nachlaufSRFull import nachlaufSRFull
from draw_conf import draw_conf

print('-------------------------------------------------------------------')
print('-                   Aufgabe 5a)                               -')
print('-------------------------------------------------------------------')


## Konstanten / Vorgaben
# (Strecken-)Lastvektor
last = np.array([0, 5])

# Anzahl Elemente
elemente = np.zeros(2, dtype="i")
elemente[0] = 10    # in x-Richtung
elemente[1] = 10    # in y-Richtung

ir = 2              # Anzahl Gausspunkte => ir^2

# Geometrie
hoehe = 44          # Hoehe
laenge = 48         # Laenge

# Materialdaten
emod = 10 ** 2      # E-Modul
nu = 0.3            # Querkontraktionszahl (Poissonzahl)
# nu = 0.4999

# Materialmatrix
C = emod/((1+nu)*(1-2*nu))*np.array([[1-nu, nu, 0],
                                     [nu, 1-nu, 0],
                                     [0, 0, (1-2*nu)/2]])


## Netzgenerator
# Netz generieren (Elementfreiheitsgradzuordnungstabelle, Geometrievektor)
edofX, x = netz(elemente[0], elemente[1], laenge, hoehe)
alphaDofsAll = elemente[0]*elemente[1]*4
alpha = np.zeros((alphaDofsAll))
edofAlpha = np.reshape(np.arange(1,alphaDofsAll+1),(4,int(alphaDofsAll/4))).T
q=np.r_[x,alpha]
new=edofAlpha+np.max(np.max(edofX))
edof=np.c_[np.array(edofX),np.array(new)]

# Geometriedatentabellen
# bzgl. Gebiet (2D Elemente)
ed = extract(np.array(edof), np.array(q))
# bzgl. Rand (1D Elemente)
idx=len(edof)-elemente[1]
edofRandRechts = edof[int(idx):, 2:6]
edRandRechts = extract(edofRandRechts, q)
# Ausgangskonfiguration darstellen
draw(ed[:,:8])


## Assembly
#   Allozieren / Speicher reservieren
ndof = q.shape[0]
xdof = x.shape[0]
K = np.zeros((ndof,ndof))
F = np.zeros((ndof))

# Elementschleife
nel = elemente[0] * elemente[1]
dataFE = {'Fe':[None]*nel,
          'indexFi':[None]*nel,
          'KeV':[None]*nel,
          'indexKi':[None]*nel,
          'indexKj':[None]*nel}

for i in np.arange(0, nel):
    Ke, Fe = elementSRFull(ed[i,:], ir, C) # Ke muss 12,12 Fe muss 4,1
    dataFE['KeV'][i] = Ke.reshape(np.size(Ke))
    dataFE['indexKi'][i] = (np.kron(np.ones(edof.shape[1]),edof[i,:]) - np.ones(edof.shape[1]**2)).astype('i')
    dataFE['indexKj'][i] = (np.kron(edof[i,:],np.ones(edof.shape[1])) - np.ones(edof.shape[1]**2)).astype('i')
    dataFE['Fe'][i] = Fe
    dataFE['indexFi'][i] = (edof[i,:] - np.ones(edof.shape[1])).astype('i')
    
#print(Ke)
#print(Fe.shape)
# sparse assembling
K = csr_matrix((np.concatenate(dataFE['KeV']),
                (np.concatenate(dataFE['indexKi']), np.concatenate(dataFE['indexKj']))))
F = csr_matrix((np.concatenate(dataFE['Fe']),
                (np.concatenate(dataFE['indexFi']), np.zeros(np.size(np.concatenate(dataFE['indexFi'])), dtype="i"))))

## Randbedingungen einbauen
# Dirichlet-Rand
RB_1 = 0
bc = RB_1 * np.ones((elemente[1] * 2 + 2, 2))
bc[:, 0] = np.arange(0, (elemente[1] * 2 + 2)).T
# Neumann-Rand
RB_2 = last

# Randelementschleife
for i in np.arange(0, elemente[1]):
    Fe = elementNeumannRand(edRandRechts[i, 0:4], ir, RB_2)
    F = assem(F, Fe, edofRandRechts[i, 0:4])

## Gleichungsloeser
d_beta = solveq(K, F, bc)

## Darstellung Gleichgewichtslage
q = q + d_beta    # Update, Konfiguration Gleichgewichtslage
ed_neu = extract(edof, q)
draw(ed_neu[:,:8])

## Nachlaufrechnung fuer Spannungsberechnung
d_neu = extract(edof, d_beta)
s = np.zeros(int(xdof/2))			
x = np.zeros((4, nel))
y = np.zeros((4, nel))
edofSpannung = (edofX[:, 1:len(edofX):2]/2).astype('i')
H = np.zeros((int(xdof / 2), int(xdof / 2)))
for i in np.arange(0, nel):
    x[0:4, i], y[0:4, i], se, He = nachlaufSRFull(ed[i,:8], ed_neu[i,:8], d_neu[i, :], ir, C)
    s = assem(s, se, edofSpannung[i, :])
    H = assem(H, He, edofSpannung[i, :])

spannungVec = np.linalg.solve(H, s)
draw_conf(spannungVec, q, edof, edofSpannung, ed)
#plt.colorbar('jet')
