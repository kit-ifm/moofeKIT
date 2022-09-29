%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 100;
setupObject.totalTime =1;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = -1;
setupObject.newton.tolerance = 1e-8;
setupObject.integrator = 'Endpoint'; 
% setupObject.integrator = 'Midpoint'; 

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
    %für linear-elastisch:
% solidViscoObject = solidClass(dofObject);
    %für linear-viskoelastisch:
solidViscoObject = solidViscoClass(dofObject);

solidViscoObject.materialObject.name = 'HookeDamage';
solidViscoObject.elementDisplacementType = 'displacement';
E = 2.26115e+03;
nu = 3e-01;
solidViscoObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));                 % Lame parameter
solidViscoObject.materialObject.mu = E/(2*(1+nu));
solidViscoObject.materialObject.rho = 10; 

solidViscoObject.materialObject.eModul0 = 250000;
solidViscoObject.materialObject.eModul1 = 50000;
solidViscoObject.materialObject.eta1 = 10000;

nel = 10;

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

% neumannObject = neumannClass(dofObject);
% neumannObject.typeOfLoad = 'deadLoad';
% neumannObject.masterObject = solidViscoObject;
% neumannObject.forceVector = [1000];
% neumannObject.shapeFunctionObject.order = solidViscoObject.shapeFunctionObject.order;
% neumannObject.shapeFunctionObject.numberOfGausspoints = 2^(solidViscoObject.dimension-1);
% neumannObject.projectionType = 'none';
% neumannObject.timeFunction = str2func('@(t) t');
% neumannObject.meshObject.edof = edofNeumann;
% neumannObject.dimension = 1;

dirichletObject = dirichletClass(dofObject);
dirichletObject.nodeList = find(solidViscoObject.meshObject.nodes(:,1)==1);
dirichletObject.nodalDof = 1;
dirichletObject.masterObject = solidViscoObject;
%eps0=0.4;                                                                  % eps0 has to be the same in displacementHookeDamage
%omega=10;
dirichletObject.timeFunction = str2func('@(t,X) (X + 0.4*cos(10*t)) ');
dirichletObject.dimension = 1;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)


%% plot hysteresis
load("hysteresis.mat");
figure;
plot(eps_array,sigma_array);
axis;
title('Hysteresis', 'fontsize',20);
xlabel('strain $\varepsilon(t)$','Interpreter','latex', 'fontsize',15);
ylabel('stress $\sigma(t)$','Interpreter','latex', 'fontsize',15);
