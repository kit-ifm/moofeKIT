import numpy as np
from element import element


def elem_loop(i, ed, edof, ir, C):
    Ke, Fe = element(ed[i, :], ir, C)
    dataFE = {}
    dataFE['KeV'] = Ke.reshape(np.size(Ke))
    dataFE['indexKi'] = (np.kron(np.ones(edof.shape[1]), edof[i, :]) -
                          np.ones(edof.shape[1]**2)).astype('i')
    dataFE['indexKj'] = (np.kron(edof[i, :], np.ones(edof.shape[1])) -
                          np.ones(edof.shape[1]**2)).astype('i')
    dataFE['Fe'] = Fe
    dataFE['indexFi'] = (edof[i, :] - np.ones(edof.shape[1])).astype('i')

    return dataFE
