import numpy as np


def netz(elementx, elementy, laenge, hoehe):
    ## netz.py - Netzgenerator fuer bilineare isoparametrische 4-Knotenelemente
    #   edof - Zuweisung der Knoten zum jeweiligen Element
    #   q    - Lage der Knoten
    #
    #
    #     (5,6) o--------------o- - -
    #           |              |
    #           |              |
    #           |     (2)      |
    #           |              |              |
    #           |              | (9,10)       |
    #     (3,4) o--------------o--------------o- - -
    #           |              |              |
    #           |              |              |
    #           |     (1)      |      (3)     |
    #           |              |              |
    #           |              |              |
    #           o--------------o--------------o- - -
    #     (1,2)              (7,8)

    nel = elementx * elementy   # Anzahl Gesamtelemente
    k = (elementy + 1) * 2      # Anzahl Knotenfreiheitsgrade in y-Richtung
    hx = laenge / elementx      # Elementlaenge in x-Richtung

    # Allozieren
    edof = np.zeros((nel, 8), dtype="i")
    ex = np.zeros((nel, 4))
    ey = np.zeros((nel, 4))

    # Zuweisung Knotenfreiheitsgrade <-> Element
    for ii in np.arange(1, elementx+1):
        edof[np.ix_(np.arange(0, elementy)+(ii-1)*elementy, [0, 1, 6, 7])] = (np.ones((4, elementy))*(np.arange(0, elementy)*2 + (ii-1)*k)).transpose() + np.ones((elementy, 1))*np.arange(1, 5)
        edof[np.ix_(np.arange(0, elementy)+(ii-1)*elementy, [2, 3, 4, 5])] = (np.ones((4, elementy))*(np.arange(0, elementy)*2 + (ii-1)*k+k)).transpose() + np.ones((elementy, 1))*np.arange(1, 5)

    # Lage der Knoten in x - Richtung
    for jj in np.arange(1, elementx+1):
        ex[np.ix_(np.arange(0, elementy)+(jj-1)*elementy,
                  np.arange(0, 4))] = np.concatenate((hx*((np.ones((elementy, 1)))*jj-1),
                                                      hx*((np.ones((elementy, 1)))*jj),
                                                      hx*((np.ones((elementy, 1)))*jj),
                                                      hx*((np.ones((elementy, 1)))*jj-1)),
                                                      axis=1)

    # Lage der Knoten in y - Richtung
    laenge_element = laenge / elementx
    Y1 = laenge_element * (hoehe / laenge)
    Y2 = laenge_element * (16 / laenge)

    for l in np.arange(1, elementx+1):
        y1 = Y1 * (l - 1)
        y2 = Y1 * l
        y3 = Y2 * (l - 1) + hoehe
        y4 = Y2 * l + hoehe

        ya = (y3 - y1) / elementy
        yb = (y4 - y2) / elementy

        for m in np.arange(0, elementy):
            ey[m+(l-1)*elementy, 0] = y1 + ya * m
            ey[m+(l-1)*elementy, 1] = y2 + yb * m
            ey[m+(l-1)*elementy, 2] = y2 + yb * (m+1)
            ey[m+(l-1)*elementy, 3] = y1 + ya * (m+1)

    ed = np.column_stack((ex[:, 0], ey[:, 0], ex[:, 1], ey[:, 1],
                          ex[:, 2], ey[:, 2], ex[:, 3], ey[:, 3]))

    # Geometrievektor erstellen
    n = np.shape(edof)
    t = edof[:, 0:n[1]]
    q = np.zeros((2 * (elementx + 1) * (elementy + 1)))
    for i in np.arange(0, n[0]):
        q[np.ix_(t[i, :]-np.ones(n[1], dtype="i"))] = ed[i, 0:n[1]]

    return edof, q
