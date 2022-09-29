%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'block';
setupObject.saveObject.saveData = false;
setupObject.totalTime = 5;
setupObject.plotObject.flag = true;
setupObject.plotObject.view = 2;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'ExplicitEuler';
% setupObject.integrator = 'Newmark';
% setupObject.newmarkGamma         = 0.5;
% setupObject.newmarkBeta          = 0.25;
switch setupObject.integrator
    case {'Endpoint','Newmark'}        
%         setupObject.totalTimeSteps       = 100;
        setupObject.totalTimeSteps       = 10;
    case {'ExplicitEuler'}
        setupObject.totalTimeSteps       = 10000;
        %courant: dt<=4.5e-4 (12000 steps)
    otherwise
        error('not implemented')
end

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
solidObject.dimension = 2;
solidObject.elementDisplacementType = 'displacement';
solidObject.shapeFunctionObject.order = 1;
solidObject.shapeFunctionObject.numberOfGausspoints = (solidObject.shapeFunctionObject.order+1)^2;
solidObject.materialObject.name = 'HookeEVZ';
% solidObject.materialObject.name = 'NeoHooke2D';
solidObject.materialObject.rho = 1;
solidObject.materialObject.E = 1e4;
solidObject.materialObject.nu = 0.2;
[solidObject.meshObject.nodes,solidObject.meshObject.edof,boundaryEdof] = meshRectangle(50,50,10,10,solidObject.shapeFunctionObject.order);
solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + 25;

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.dimension = 2;
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find((solidObject.meshObject.nodes(:,1)==0) & (abs(solidObject.meshObject.nodes(:,2))<=10));
dirichletObject.nodalDof = 1:2;
dirichletObject.timeFunction = str2func('@(t) 0');

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.dimension = 2;
neumannObject.typeOfLoad = 'deadLoad';
neumannObject.masterObject = solidObject;
neumannObject.forceVector = [-2000;0];
neumannObject.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = neumannObject.shapeFunctionObject.order+1;
neumannObject.projectionType = 'none';
neumannObject.timeFunction = @(t) sin(pi*t).*(t <= 1);
neumannObject.meshObject.edof = boundaryEdof.SX2;

%% solver
dofObject = runNewton(setupObject,dofObject);
% plot(solidObject)

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end,2)-solidObject.qR(end,2);
fprintf('\nDisplacement: %4.3f\n',yDisplacementUpperRightNode)