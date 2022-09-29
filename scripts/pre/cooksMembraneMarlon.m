%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembraneMarlon';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 4;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = -1;
setupObject.plotObject.postPlotType = 'stress';
setupObject.newton.tolerance = 1e-4;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.materialObject.name = 'Hooke';
solidObject.elementDisplacementType = 'displacement';

E = 2.26115e+03;
nu = 3e-01;
solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
solidObject.materialObject.mu = E/(2*(1+nu));
solidObject.materialObject.rho = 0;
nel = 4;
[solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,nel,48,44,16,10);
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = 8;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1:3;
dirichletObject.masterObject = solidObject;
dirichletObject.timeFunction = str2func('@(t) 0');

neumannObject = neumannClass(dofObject);
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidObject;
neumannObject.forceVector = [0;50;0];
% neumannObject.forceVector = [-50;0;0];
neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject.projectionType = 'none';
neumannObject.timeFunction = str2func('@(t) t');
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
plot(solidObject,setupObject)
