import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from extract import extract

def draw_conf(Spann,q,edof,edofSpannung,ed):
    
    ## Konturplot Spannung

    ed_neu = extract(edof, q)
    ed_Spann = extract(edofSpannung, Spann)
    
    X = np.zeros((4,np.shape(ed)[0]))
    Y = np.zeros((4,np.shape(ed)[0]))
    PK = np.zeros((4,np.shape(ed)[0]))
    for i in np.arange(0,np.shape(ed)[0]):
        X[0:4,i] = ed_neu[np.ix_([i], [0, 2, 4, 6])]
        Y[0:4,i] = ed_neu[np.ix_([i], [1, 3, 5, 7])]
        PK[0:4,i] = ed_Spann[np.ix_([i], [0, 1, 2, 3])]
    
    # ptc.Polygon(X, Y, PK)
    # caxis(mcat([0, max(Spann)]))
    # shading interp
    # set(gcf, mstring('color'), mcat([1, 1, 1]))