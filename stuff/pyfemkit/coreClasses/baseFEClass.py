import numpy as np
from commonFunctions.checkIntegerMatrix import checkIntegerMatrix

class baseFEClass():
    def __init__(self):
        self.nodes = np.array([],dtype=float)
        self.edof = np.array([],dtype=int)             # local node based edof (element degrees of freedom table)
        self.globalFullEdof = np.array([],dtype=int)     # global dof based edof
        self.globalNodesDof = np.array([],dtype=int)    # global node based dof (degrees of freedom)
        self.dimension = 0
#set methods
    def setedof(self,input): #set. nicht möglich
        checkIntegerMatrix(input)
        self.edof=int(input)

    def setglobalFullEdof(self,input):
        checkIntegerMatrix(input)
        self.globalFullEdof=int(input)

    def setglobalNodesDof(self,input):
        checkIntegerMatrix(input)
        self.globalNodesDof = int(input)
        
    def setdimension(self,input):
        checkIntegerMatrix(input)
        if (input>3 and input <1):
            self.dimension = input
