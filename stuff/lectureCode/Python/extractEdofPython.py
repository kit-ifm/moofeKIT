import numpy as np


def extractEdofPython(edof, q):
    # Extrahieren der Elementdatentabelle aus der
    # Elementfreiheitsgradzuordnungsdatentabelle und dem Geometrievektor q
    #
    # ed:   Elementdatentabelle
    # edof: Elementfreiheitsgradzuordnungsdatentabelle
    # q:    Geometrievektor

    ed = q[edof.T].T

    return ed
