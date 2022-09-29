import numpy as np
from scipy.sparse import csr_matrix
from commonFunctions.lagrangeShapeFunctions import lagrangeShapeFunctions

class dofClass(object):
    def __init__(self):
        self.qN = np.array([],dtype=float)
        self.qN1 = np.array([],dtype=float)
        self.vN = np.array([],dtype=float)
        self.vN1 = np.array([],dtype=float)
        self.listContinuumObjects = []
        self.totalNumberOfDofs = 0
        # self.dirichletDof
        self.numberOfContinuumObjects
        # self.solveDof
        self.callMassMatrixElement

    @property
    def numberOfContinuumObjects(self):
        return len(self.listContinuumObjects)
        
    # @property
    def dirichletDof(self):
        dirichletDofArray = []
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            if hasattr(continuumObject,'globalDirichletList'):  
                dirichletDofArray = np.unique(np.concatenate((dirichletDofArray,continuumObject.globalDirichletList()), axis = None))
        return dirichletDofArray
    
    # @property
    def solveDof(self):
        solveDofArray = list(range(self.totalNumberOfDofs))
        if len(self.dirichletDof())!=0:
            for index1 in range(len(self.dirichletDof())):
                solveDofArray.remove(self.dirichletDof()[index1])
            return solveDofArray

    def initialize(self,fieldNameCell):
        for index1 in range(len(fieldNameCell)):
            fieldName = fieldNameCell[index1]
            assert isinstance(fieldName, str),'input need to be of type string'  
            setattr(self,fieldName,np.zeros(self.totalNumberOfDofs))

    def updateGlobalField(self,fieldNameCell):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            for index2 in range(len(fieldNameCell)):
                fieldName = fieldNameCell[index2]
                assert isinstance(fieldName, str),'input need to be of type string'               
                if hasattr(continuumObject,fieldName):
                    getattr(self,fieldName)[continuumObject.globalNodesDof] = getattr(continuumObject,fieldName)
            ## Test
            # self.vN1 = np.array([1,1,1.0,1,1,1])
            # fieldName2 = 'vN1'
            # ind1 = np.array([1,3,5],int)
            # ind2 = np.array([0,1,2],int)
                # if hasattr(self,fieldName2):
                    # getattr(self,fieldName)[ind1] = getattr(self,fieldName2)[ind2]

    def updateTimeDependentField(self,fieldNameCell,time,flag):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            for index2 in range(len(fieldNameCell)):
                fieldName = fieldNameCell[index2]
                assert isinstance(fieldName, str),'input need to be of type string'               
                if hasattr(continuumObject,fieldName):
                    if continuumObject.flagNewton[flag]:
                        continuumObject.time = time
                        getattr(self,fieldName)[continuumObject.globalDirichletList()] = getattr(continuumObject,fieldName)()
                        
    def updateContinuumField(self,fieldNameCell,flag):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            for index2 in range(len(fieldNameCell)):
                fieldName = fieldNameCell[index2]
                assert (isinstance(fieldName, str) and isinstance(flag, str)),'input need to be of type string'
                field = getattr(self,fieldName)
                if (hasattr(continuumObject,fieldName) and flag in continuumObject.flagNewton):
                    if continuumObject.flagNewton[flag]:
                        setattr(continuumObject,fieldName,field[continuumObject.globalNodesDof])

    @property
    def callMassMatrixElement(self):
        # massmassDataFEData = []
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            if (hasattr(continuumObject,'massMatrix')):  
                # TODO
                # massData = [massData, continuumObject.massMatrix()]
                massDataFE = continuumObject.massMatrix()
            # assembly
            return csr_matrix((np.concatenate(massDataFE['MeVector']),(np.concatenate(massDataFE['indexMi']), np.concatenate(massDataFE['indexMj'])))) 
    
    def initializeGlobalDofs(self,flag):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            if flag in continuumObject.flagNewton:
                if continuumObject.flagNewton[flag]:
                    numberOfDofPerNode = continuumObject.nodes.shape[1]
                    edof = continuumObject.edof
                    numberOfNodes = continuumObject.nodes.shape[0]
                    numberOfElements = edof.shape[0]
                    globalFullEdof = np.zeros((numberOfElements,edof.shape[1]*numberOfDofPerNode))
                    for j in range(numberOfDofPerNode):
                        globalFullEdof[:,j::numberOfDofPerNode] = self.totalNumberOfDofs + (edof+1)*numberOfDofPerNode - (numberOfDofPerNode-(j+1))
                    globalFullEdof = np.array(globalFullEdof - 1,dtype=int)
                    # globalFullEdof(:,j:numberOfDofPerNode:end) = obj.totalNumberOfDofs + edof*numberOfDofPerNode-(numberOfDofPerNode-j);
                    globalNodesDof = np.zeros((numberOfNodes,numberOfDofPerNode))
                    globalNodesDof[:,0:numberOfDofPerNode+1] = self.totalNumberOfDofs + np.kron(np.ones(numberOfNodes).reshape(-1,1),range(1,numberOfDofPerNode+1)) + np.kron(np.arange(0,numberOfDofPerNode*(numberOfNodes),numberOfDofPerNode).reshape([-1,1]),np.ones(numberOfDofPerNode))
                    #
                    globalNodesDof = np.array(globalNodesDof - 1,dtype=int)
                    #
                    continuumObject.globalFullEdof = globalFullEdof
                    continuumObject.globalNodesDof = globalNodesDof
                    self.totalNumberOfDofs = self.totalNumberOfDofs + numberOfNodes*numberOfDofPerNode

    def initializeQV(self,flag):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            if flag in continuumObject.flagNewton:
                if continuumObject.flagNewton[flag]:
                    if continuumObject.qN.size == 0:
                        continuumObject.qN = continuumObject.qR
                    if continuumObject.qN1.size == 0:
                        continuumObject.qN1 = continuumObject.qR
                    if continuumObject.vN.size == 0:
                        continuumObject.vN = np.zeros(continuumObject.nodes.shape)
                        
    def initializeShapeFunctions(self,flag):
        for index1 in range(self.numberOfContinuumObjects):
            continuumObject = self.listContinuumObjects[index1]
            if flag in continuumObject.flagNewton:
                if continuumObject.flagNewton[flag]:
                    continuumObject.shapeFunctions = lagrangeShapeFunctions(continuumObject.edof.shape[1], continuumObject.numberOfGausspoints, continuumObject.dimension)
