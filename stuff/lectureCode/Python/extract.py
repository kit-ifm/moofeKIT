import numpy as np


def extract(edof, q):
    # Extrahieren der Elementdatentabelle aus der
    # Elementfreiheitsgradzuordnungsdatentabelle und dem Geometrievektor q
    #
    # ed:   Elementdatentabelle
    # edof: Elementfreiheitsgradzuordnungsdatentabelle
    # q:    Geometrievektor

    ed = q[(edof - np.ones(edof.shape, dtype="i")).T].T

    return ed
