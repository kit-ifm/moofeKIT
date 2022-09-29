%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 40;
setupObject.totalTime = 10;
setupObject.plotObject.flag = false;
setupObject.plotObject.stress.component = -1;
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Endpoint'; 

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
    %für linear-elastisch:
solidViscoObject = solidClass(dofObject);
solidViscoObject.materialObject.name = 'Hooke';
solidViscoObject.elementDisplacementType = 'displacement';
E = 2.26115e+03;
nu = 3e-01;
solidViscoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));                 % Lame parameter
solidViscoObject.materialObject.mu = E/(2*(1+nu));
solidViscoObject.materialObject.rho = 0; 

nel = 2;
[solidViscoObject.meshObject.nodes,solidViscoObject.meshObject.edof,edofNeumann] = netzdehnstab(nel,1);

solidViscoObject.dimension = 1;
solidViscoObject.shapeFunctionObject.order = 1;
solidViscoObject.shapeFunctionObject.numberOfGausspoints = 2;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1;
dirichletObject.masterObject = solidViscoObject;
dirichletObject.timeFunction = str2func('@(t,X) 0*X');
dirichletObject.dimension = 1;

neumannObject = neumannClass(dofObject);
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidViscoObject;
neumannObject.forceVector = [1000];
neumannObject.shapeFunctionObject.order = solidViscoObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidViscoObject.dimension-1);
neumannObject.projectionType = 'none';
neumannObject.timeFunction = str2func('@(t) t');
neumannObject.meshObject.edof = find(solidViscoObject.meshObject.nodes(:,1)==0.5);
neumannObject.dimension = 1;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==1);
dirichletObject.nodalDof = 1;
dirichletObject.masterObject = solidViscoObject;
dirichletObject.timeFunction = str2func('@(t,X) 0*X');
dirichletObject.dimension = 1;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)