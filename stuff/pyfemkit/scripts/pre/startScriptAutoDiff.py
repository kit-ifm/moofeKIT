# libraries
import numpy as np
import matplotlib.pyplot as plt
from ttictoc import tic,toc

# import classes & functions
from coreClasses.setupClass import setupClass
from coreClasses.dofClass import dofClass
from continuumClasses.solidClass.solidClass import solidClass
from continuumClasses.dirichletClass.dirichletClass import dirichletClass
from meshes.meshGenerators.meshGeneratorCube import meshGeneratorCube
from scripts.solver.runNewton import runNewton

# setupObject (mandatory)
setupObject = setupClass()
setupObject.fileName = 'patchTest'
setupObject.totalTimeSteps = 1
setupObject.totalTime = 1.0
setupObject.plotFlag = False
# dofObject (mandatory)
dofObject = dofClass()
# continuumObjects
part1 = solidClass()
part1.solidClassfunc(dofObject)

# part1.materialName = 'Hooke'
# part1.materialName = 'HookeAutoDiff'
# part1.materialName = 'SaintVenant'
part1.materialName = 'SaintVenantAutoDiff'
part1.materialData['Lambda'] = 27.7777
part1.materialData['mu'] = 41.6666
part1.materialData['rho'] = 1
part1.nodes,part1.edof,bounEDOFs = meshGeneratorCube(1,1,1,5,5,5,1,False)
part1.edof = part1.edof - 1
part1.nodes = part1.nodes + 0.5
part1.dimension = 3
part1.orderShapeFunctions = 1
part1.numberOfGausspoints = 8
# dirichletObjects
boundary1 = dirichletClass()
boundary1.dirichletClassfunc(dofObject)
boundary1.nodeList = np.array(np.where(part1.nodes[:,2] == 1),dtype=int)
boundary1.nodalDof = 3
boundary1.masterObject = part1
# solver
# tic()
dofObject = runNewton(setupObject,dofObject)
# print(toc())
part1.plot('pyvista')
# print(dofObject.listContinuumObjects[-1].__dict__)
