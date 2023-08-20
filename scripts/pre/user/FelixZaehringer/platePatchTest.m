%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 2;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type','Cauchy', 'component', 13);
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof, edofBoundary] = meshPatchTestDistorted2D(0.24, 0.12, 1, false);
plateObject.materialObject.name = 'HookeBD';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e6;
plateObject.materialObject.nu = 0.25;
plateObject.h = 0.001;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = 1;
plateObject.shapeFunctionObject.numberOfGausspoints = (1+1)^2;
plateObject.elementDisplacementType = 'displacement';
% plateObject.mixedFEObject.condensation = true;

% Dirichlet boundary
nodes = plateObject.meshObject.nodes;
numberOfNodes = size(nodes, 1);

for ii=1:numberOfNodes
    x = num2str(nodes(ii, 1));
    y = num2str(nodes(ii, 2));

    boundary2 = dirichletClass(dofObject);
    boundary2.nodeList = ii;
    boundary2.nodalDof = 1;
    boundary2.timeFunction = str2func(['@(t, nodalCoordinates) t*1e-3*(', x, '^2+',y,'^2+',x,'*',y,')/2']);
    boundary2.masterObject = plateObject;

    boundary3 = dirichletClass(dofObject);
    boundary3.nodeList = ii;
    boundary3.nodalDof = 2;
    boundary3.timeFunction = str2func(['@(t, nodalCoordinates) t*1e-3*(',x,'/2+',y,')']);
    boundary3.masterObject = plateObject;

    boundary4 = dirichletClass(dofObject);
    boundary4.nodeList = ii;
    boundary4.nodalDof = 3;
    boundary4.timeFunction = str2func(['@(t, nodalCoordinates) t*1e-3*(-',x,'-',y,'/2)']);
    boundary4.masterObject = plateObject;
end


%% solver
dofObject = runNewton(setupObject,dofObject);
%plot(plateObject, setupObject)
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);