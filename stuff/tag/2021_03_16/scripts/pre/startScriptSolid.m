run('../prepareWorkspace');

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.fileName = 'solidDiscreteGradient';
setupObject.totalTimeSteps = 4;
setupObject.totalTime = 1;
setupObject.plotFlag = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
part2 = solidClass(dofObject);
[part2.nodes,part2.edof] = meshGeneratorCube(1,1,1,3,3,3,1,false);
part2.nodes = part2.nodes + 2.5;
part2.elementDisplacementType = 'displacementSC';
part2.materialName = 'MooneyRivlin';
setupObject.integrator = 'DiscreteGradient';
part2.materialData.rho = 1;
part2.materialData.a = 1;
part2.materialData.b = 1;
part2.materialData.c = 1;
part2.materialData.d = 2*(part2.materialData.a + 2*part2.materialData.b);
part2.dimension = 3;
part2.orderShapeFunctions = 1;
part2.numberOfGausspoints = 8;
 
boundaryK2 = dirichletClass(dofObject);
boundaryK2.nodeList = find(part2.nodes(:,3) == 3);
boundaryK2.nodalDof = 3;
boundaryK2.masterObject = part2;
boundaryK2.timeFunction = str2func('@(t,Z) (Z - 0.25*t).*(t >= 0)');

%% solver
dofObject = runNewton(setupObject,dofObject);