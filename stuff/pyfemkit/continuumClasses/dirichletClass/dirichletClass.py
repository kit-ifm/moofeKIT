import numpy as np
import inspect
import sys
from commonFunctions.checkIntegerMatrix import checkIntegerMatrix

class dirichletClass(object):
    def __init__(self):
        self.flagNewton = {'initializeGlobalDofs': False,
                            'initializeQV': False,
                            'initializeShapeFunctions': False,
                            'updateTimeDependentField': True,
                            'updateContinuumFieldPreNewtonLoop': False,
                            'callElements': False,
                            'updateContinuumFieldNewtonLoop': False}
        #properties        
        self.masterObject
        self._nodeList = np.array([],dtype=int)
        self._nodalDof = 1   # 1,2,3  
        self._time = 0
        #dependent=true
        # self._globalDirichletList = np.array([],dtype=float) 
        # self.qN1 # = np.array([],dtype=float) 

    def dirichletClassfunc(self,dofObject):
        dofObject.listContinuumObjects.append(self) 
    
    def masterObject(self,input):
        if (inspect.isclass(input,'solidClass') or inspect.isclass(input,'solidThermoClass')):
            self.masterObject = input
        else:
            sys.exit('masterObject must be compatibel to dirichlet class!')
        
    def timeFunction(self,time,qR):
        # linear function
        # t = qR-time
        t = qR-0.5*time
        return t

    @property
    def nodeList(self):
        return self._nodeList
    @nodeList.setter
    def nodeList(self, input):
        # checkIntegerMatrix(input)
        self._nodeList = np.int32(input.reshape(-1,1)) 

    @property
    def nodalDof(self):
        return self._nodalDof
    @nodalDof.setter
    def nodalDof(self, input):
        # checkIntegerMatrix(input)solveDof
        self._nodalDof = input
            
    @property
    def time(self):
        return self._time
    @time.setter
    def time(self, input):
        assert isinstance(input, float),'time need to be class double!'
        self._time = input

    # @property
    def globalDirichletList(self):
        if not(self.masterObject.globalNodesDof.size == 0):
            lists = self.masterObject.globalNodesDof[self.nodeList,self.nodalDof-1]
            return lists.reshape(-1,1)
        else:
            return np.array([],dtype=int)
        
    # @property
    def qN1(self):
        if np.size(self.nodalDof)==1:
            return self.timeFunction(self.time,self.masterObject.qR[self.nodeList,self.nodalDof-1])
        else:
            return self.masterObject.qR[self.nodeList,self.nodalDof-1] + self.timeFunction[self.time]
