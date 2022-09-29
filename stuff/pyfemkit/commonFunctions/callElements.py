import numpy as np
from continuumClasses.solidClass.displacementHookeEndpoint import displacementHookeEndpoint
from continuumClasses.solidClass.displacementHookeAutoDiffEndpoint import displacementHookeAutoDiffEndpoint
from continuumClasses.solidClass.displacementParallelHookeAutoDiffEndpoint import displacementParallelHookeAutoDiffEndpoint
from continuumClasses.solidClass.displacementSaintVenantEndpoint import displacementSaintVenantEndpoint
from continuumClasses.solidClass.displacementSaintVenantAutoDiffEndpoint import displacementSaintVenantAutoDiffEndpoint
from continuumClasses.solidClass.displacementSCMooneyRivlinDiscreteGradient import displacementSCMooneyRivlinDiscreteGradient
from continuumClasses.solidClass.displacementSCMooneyRivlinAutoDiffDiscreteGradient import displacementSCMooneyRivlinAutoDiffDiscreteGradient
from continuumClasses.solidClass.displacementParallelHookeEndpoint import displacementParallelHookeEndpoint

def callElements(self,setupObject,flag):
    dataFE = []
    for index1 in range(self.numberOfContinuumObjects):
        continuumObject = self.listContinuumObjects[index1]
        if continuumObject.flagNewton[flag]:
            elementName = continuumObject.elementDisplacementType + continuumObject.materialName + setupObject.integrator
            dataFEindex = eval(elementName)(continuumObject,setupObject)
            # FIXME
            # dataFE = [dataFE, dataFEindex]
            dataFE = dataFEindex
    return dataFE