import numpy as np
import sys

class setupClass(object):
    def __init__(self):
        self._fileName = 'save'
        self._newtonTolerance = 1e-8
        self._newtonMaximumSteps = 50
        self._newtonMaximumResiduumNorm = 1e12
        self._totalTime = 1.0
        self._totalTimeSteps = 1
        self._integrator = 'Endpoint'
        self.timeStepSize
        self.factorIntegrator      
        # TODO
        self.meshType = 'meshGenerator' #abaqusInput
        self.meshName = 'Square'
        #self.timeFunction #= @(t) obj.timeStepSize
        self.plotFlag = False
        self.makeMovie = False

    @property
    def fileName(self):
        return self._fileName
    @fileName.setter
    def fileName(self, input):
        assert isinstance(input, str),'input need to be of type string'
        self._fileName = input
        
    @property
    def newtonTolerance(self):
        return self._newtonTolerance
    @newtonTolerance.setter
    def newtonTolerance(self, input):
        assert (isinstance(input, float) and 0<input<1),'input need to be of type float and inbetween 0 and 1'
        self._newtonTolerance = input
        
    @property
    def newtonMaximumSteps(self):
        return self._newtonMaximumSteps
    @newtonMaximumSteps.setter
    def newtonMaximumSteps(self, input):
        assert (isinstance(input, int) and input>0),'input need to be of type integer and greater 0'
        self._newtonMaximumSteps = input

    @property
    def newtonMaximumResiduumNorm(self):
        return self._newtonMaximumResiduumNorm
    @newtonMaximumResiduumNorm.setter
    def newtonMaximumResiduumNorm(self, input):
        assert (isinstance(input, float) and input>0),'input need to be of type float and greater 0'
        self._newtonMaximumResiduumNorm = input
                    
    @property
    def totalTime(self):
        return self._totalTime
    @totalTime.setter
    def totalTime(self, input):
        assert (isinstance(input, float) and input>0),'input need to be of type float and greater 0'
        self._totalTime = input

    @property
    def totalTimeSteps(self):
        return self._totalTimeSteps
    @totalTimeSteps.setter
    def totalTimeSteps(self, input):
        assert (isinstance(input, int) and input>0),'input need to be of type integer and greater 0'
        self._totalTimeSteps = input
               
    @property
    def integrator(self):
        return self._integrator
    @integrator.setter
    def integrator(self, input):
        assert isinstance(input, str),'input need to be of type string'
        self._integrator = input

# get methods
    @property
    def timeStepSize(self):
        return self._totalTime/self._totalTimeSteps

    @property
    def factorIntegrator(self):
        if self.integrator == 'Endpoint':
            return [1, 0] #np.array([1, 0],dtype=int)
        elif (self.integrator == 'Endpoint' or self.integrator == 'DiscreteGradient'):
            return [2, 2] #np.array([2, 2],dtype=int)
        else:
            sys.exit('integrator is not available')

