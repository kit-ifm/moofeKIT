%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'cooksMembrane';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 4;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = 11;
setupObject.newton.tolerance = 1e-4;
% setupObject.integrator = 'Endpoint';
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';
% 
dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
% solidObject.materialObject = materialPackage.saintVenantClass;
solidObject.materialObject = materialPackage.neoHookeClass;
% solidObject.materialObject = materialPackage.neoHookeSplitClass; solidObject.materialObject.kappa = 1;
% solidObject.materialObject = materialPackage.logStrainClass; solidObject.materialObject.kappa = 1; solidObject.materialObject.specDecompAlgo = 'exactEigenVec'; solidObject.materialObject.perturbThreshold = 1e5*eps;
% solidObject.materialObject.E = 2.26115e+03;
% solidObject.materialObject.nu = 4.95469e-01;
solidObject.materialObject.materialDatabase('academic');
% solidObject.materialObject = materialPackage.mooneyRivlinClass; 
% shear   = 997.5;mu      = shear;solidObject.materialObject.c1      = 5/6*mu;solidObject.materialObject.c2      = 1/6*mu;solidObject.materialObject.c       = 0;
% solidObject.materialObject = materialPackage.ogdenClass;
% solidObject.materialObject.mue = 997.5; solidObject.materialObject.alpha = 1; solidObject.materialObject.beta = 1; solidObject.materialObject.kappa = 1; solidObject.materialObject.perturbThreshold = 1e5*eps; solidObject.materialObject.specDecompAlgo = 'exactEigenVec';%'perturbation' or 'exacteigenvec'
solidObject.materialObject.rho = 10;
nel = 2;
[solidObject.meshObject.nodes,solidObject.meshObject.edof,edofNeumann] = trilinearCooksMembrane(nel,nel,nel,48,44,16,10);
% [solidObject.nodes,solidObject.edof,edofNeumann] = linearTetrahedralCooksMembrane(nel,nel,nel,48,44,16,10);
% FIXME: negative determinants for quadratic 27-node Element
% [solidObject.nodes,solidObject.edof,edofNeumann] = triquadraticCooksMembrane(nel,nel,nel,48,44,16,10);
% FIXME: negative determinants for quadratic 10-node Element
% [solidObject.nodes,solidObject.edof,edofNeumann] = quadraticTetrahedralCooksMembrane(nel,nel,nel,48,44,16,10);
% solidObject.materialObject.name = 'NeoHooke';
solidObject.materialObject.name = 'Hooke';
% solidObject.materialObject.name = 'SaintVenant';
% solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
% solidObject.materialObject.name = 'HyperelasticPF'; % 'Hooke';
% solidObject.elementDisplacementType = 'easSC'; %'displacement'; % 'displacementSC'; % 'eas';
solidObject.elementDisplacementType = 'displacement'; % 'displacementSC'; % 'eas';
% solidObject.materialObject.rho = 10;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 2;
solidObject.shapeFunctionObject.numberOfGausspoints = 27;
solidObject.mixedFEObject.condensation = false;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:,1)==0);
dirichletObject.nodalDof = 1:3;
dirichletObject.masterObject = solidObject;
dirichletObject.timeFunction = str2func('@(t) 0');

neumannObject = neumannClass(dofObject);
neumannObject.dimension = 3;
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidObject;
neumannObject.forceVector = [0;50;0];
neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject.projectionType = 'none';
neumannObject.timeFunction = str2func('@(t) t');
neumannObject.meshObject.edof = edofNeumann;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)