import numpy as np
import matplotlib.pyplot as plt
# plot 
import pyvista as pv
import vtk
#from coreClasses.baseFEClass import baseFEClass
from coreClasses.solidSuperClass.solidSuperClass import solidSuperClass
# from structClass import structClass
from commonFunctions.checkIntegerMatrix import checkIntegerMatrix
from coreClasses.solidSuperClass.massMatrixElement import massMatrixElement

class solidClass(solidSuperClass):#object):
    def __init__(self):
        self.filename = 'save'
        self.flagNewton = {'initializeGlobalDofs': True,
                            'initializeQV': True,
                            'initializeShapeFunctions': True,
                            'updateTimeDependentField': False,
                            'updateContinuumFieldPreNewtonLoop': True,
                            'callElements': True,
                            'updateContinuumFieldNewtonLoop': True}
        self.additionalFields = 0
        self.elementDisplacementType = 'displacement'
        self.plot
        solidSuperClass.__init__(self)

    def solidClassfunc(self,dofObject):
        dofObject.listContinuumObjects.append(self) 
        
    def plot(self,input):
        if input == 'pyvista':
            cellsVector = (self.edof.shape[1]*np.ones(self.edof.shape[0])).astype(int)
            cells = self.edof.copy()
            # col1 = 5
            # col2 = 2
            # cells.T[[col1, col2]] = cells.T[[col2, col1]]
            # col1 = 4
            # col2 = 3
            # cells.T[[col1, col2]] = cells.T[[col2, col1]]
            # col1 = 2
            # col2 = 5
            # cells.T[[col1, col2]] = cells.T[[col2, col1]]
            # col1 = 3
            # col2 = 4
            # cells.T[[col1, col2]] = cells.T[[col2, col1]]
            cells = np.append(cellsVector.reshape(-1,1),cells,axis=1)
            cells = cells.ravel()
            celltypes = np.zeros(self.edof.shape[0], dtype=np.uint8)
            celltypes[:] = vtk.VTK_HEXAHEDRON
            grid = pv.UnstructuredGrid(cells, celltypes, self.qN1)
            _ = grid.plot(show_edges=True)
        elif input == 'matlab': 
            import matlab.engine
            eng = matlab.engine.start_matlab()
            SX = np.array([[1, 2, 3, 4],
                           [5, 8, 7, 6],
                           [1, 5, 6, 2],
                           [2, 6, 7, 3],
                           [3, 7, 8, 4],
                           [4, 8, 5, 1]],dtype=int)
            SX = SX - np.ones(SX.shape, dtype="i")
            Faces = np.concatenate([self.edof[:,SX[0,:]], self.edof[:,SX[1,:]], self.edof[:,SX[2,:]], self.edof[:,SX[3,:]], self.edof[:,SX[4,:]], self.edof[:,SX[5,:]]]) + 1
            FaceColor = np.array([0,1,0])
            FaceAlpha = 0.5
            eng.patch('Vertices',matlab.double(self.qN1.tolist()),'Faces',matlab.uint32(Faces.tolist()),'FaceColor',matlab.double(FaceColor.tolist()),'FaceAlpha',FaceAlpha)
            eng.eval('view(3)')
            eng.eval('axis on')