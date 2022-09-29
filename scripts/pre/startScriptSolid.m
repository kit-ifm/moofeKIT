%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'solid';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
[solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,1,1,1,1,1,1,false);%
solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
solidObject.elementDisplacementType = 'mixedSC';
solidObject.materialObject.name = 'MooneyRivlin';
setupObject.integrator = 'DiscreteGradient';
solidObject.materialObject.rho = 1;
solidObject.materialObject.a = 1;
solidObject.materialObject.b = 1;
solidObject.materialObject.c = 1;
solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;
solidObject.mixedFEObject.typeShapeFunctionData = 1;
solidObject.mixedFEObject.condensation = false;
solidObject.numericalTangentObject.computeNumericalTangent = false;
solidObject.numericalTangentObject.showDifferences = true;
solidObject.numericalTangentObject.type = 'complex';

dirichletObject1 = dirichletClass(dofObject);
dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:,3) == 1);
dirichletObject1.nodalDof = 3;
dirichletObject1.masterObject = solidObject;
dirichletObject1.timeFunction = str2func('@(t,Z) (Z - 0.5).*(t >= 1) + (Z - 0.5*t).*(t >= 0).*(t < 1)');

dirichletObject2 = dirichletClass(dofObject);
dirichletObject2.nodeList = find(solidObject.meshObject.nodes(:,3) == 0);
dirichletObject2.nodalDof = 3;
dirichletObject2.masterObject = solidObject;

dirichletObject3 = dirichletClass(dofObject);
dirichletObject3.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
dirichletObject3.nodalDof = 1;
dirichletObject3.masterObject = solidObject;

dirichletObject4 = dirichletClass(dofObject);
dirichletObject4.nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
dirichletObject4.nodalDof = 2;
dirichletObject4.masterObject = solidObject;

%% solver
dofObject = runNewton(setupObject,dofObject);