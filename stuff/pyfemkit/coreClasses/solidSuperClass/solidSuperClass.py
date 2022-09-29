import numpy as np
from coreClasses.solidSuperClass.massMatrixElement import massMatrixElement
from coreClasses.baseFEClass import baseFEClass

class solidSuperClass(baseFEClass):
    flagNewton = None
    additionalField = None # constant and abstract außerhalb __init__ Umgebung
    def __init__(self):
        #abstract and constant --> zu klären. Als constant-Element nicht möglich zwecks Vererbung
        #self.flagNewton
        #self.additionalField

        #verberbung
        baseFEClass.__init__(self) 

        #mesh
        self.orderShapeFunctions = 1 
        self.numberOfGausspoints = 1
        self.shapeFunctions =[]

        #material
        self.materialData = {}
        self.materialName = None
        # self.materialData = None

        #solver
        self.qN = np.array([],dtype=float)
        self.qN1 = np.array([],dtype=float)
        self.vN = np.array([],dtype=float)
        self.vN1 = np.array([],dtype=float)

        #dependent
        self.qR
        self.massMatrix

#get methods
    @property
    def qR(self): #get. nicht möglich
        return self.nodes
    
    #@property
    def massMatrix(self):
        return massMatrixElement(self)

        
