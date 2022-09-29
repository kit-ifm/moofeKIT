## Importieren benoetigter Bibliotheken und Funktionen
# numPy & SciPy
import numpy as np
from scipy.sparse import csr_matrix

# Funktionen
from netz import netz
from extract import extract
from draw import draw
#from element import element
from elementSRReduced import elementSRReduced
from assem import assem
from elementNeumannRand import elementNeumannRand
from solveq import solveq
#from nachlauf import nachlauf
from nachlaufSRReduced import nachlaufSRReduced
from draw_conf import draw_conf

print('-------------------------------------------------------------------')
print('-                   Aufgabe 5b)                               -')
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
edof, q = netz(elemente[0], elemente[1], laenge, hoehe)

# Geometriedatentabellen
# bzgl. Gebiet (2D Elemente)
ed = extract(edof, q)
# bzgl. Rand (1D Elemente)
edofRandRechts = edof[-elemente[1]:, 2:6]
edRandRechts = extract(edofRandRechts, q)

# Ausgangskonfiguration darstellen
draw(ed)


## Assembly
#   Allozieren / Speicher reservieren
ndof = q.shape[0]
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
    Ke, Fe = elementSRReduced(ed[i, :], ir, C)
    dataFE['KeV'][i] = Ke.reshape(np.size(Ke))
    dataFE['indexKi'][i] = (np.kron(np.ones(edof.shape[1]),edof[i,:]) - np.ones(edof.shape[1]**2)).astype('i')
    dataFE['indexKj'][i] = (np.kron(edof[i,:],np.ones(edof.shape[1])) - np.ones(edof.shape[1]**2)).astype('i')
    dataFE['Fe'][i] = Fe
    dataFE['indexFi'][i] = (edof[i,:] - np.ones(edof.shape[1])).astype('i')
    

# sparse assembling
K = csr_matrix((np.concatenate(dataFE['KeV']),
                (np.concatenate(dataFE['indexKi']), np.concatenate(dataFE['indexKj']))))
F = csr_matrix((np.concatenate(dataFE['Fe']),
                (np.concatenate(dataFE['indexFi']), np.zeros(np.size(np.concatenate(dataFE['indexFi'])), dtype="i"))))

## Randbedingungen einbauen
# Dirichlet-Rand
RB_1 = 0
bc = RB_1 * np.ones((elemente[1] * 2 + 2, 2))
bc[:, 0] = np.arange(1, (elemente[1] * 2 + 3)).T
# Neumann-Rand
RB_2 = last

# Randelementschleife
for i in np.arange(0, elemente[0]):
    Fe = elementNeumannRand(edRandRechts[i, :], ir, RB_2)
    F = assem(F, Fe, edofRandRechts[i, :])

## Gleichungsloeser
d = solveq(K, F, bc)

## Darstellung Gleichgewichtslage
q = q + d    # Update, Konfiguration Gleichgewichtslage
ed_neu = extract(edof, q)
draw(ed_neu)

## Nachlaufrechnung fuer Spannungsberechnung
d_neu = extract(edof, d)
s = np.zeros(int(ndof / 2))						
x = np.zeros((4, nel))
y = np.zeros((4, nel))
edofSpannung = (edof[:, 1::2] / 2).astype('i')
H = np.zeros((int(ndof / 2), int(ndof / 2)))
for i in np.arange(0, nel):
    x[0:4, i], y[0:4, i], se, He = nachlaufSRReduced(ed[i, :], ed_neu[i, :],	d_neu[i, :], ir, C)
    s = assem(s, se, edofSpannung[i, :])
    H = assem(H, He, edofSpannung[i, :])

spannungVec = np.linalg.solve(H, s)
draw_conf(spannungVec, q, edof, edofSpannung, ed)
#plt.colorbar('jet')
