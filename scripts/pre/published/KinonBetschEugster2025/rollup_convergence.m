%% Static GE beam with point torque at free end
clearvars; close all;

%% Parameters
lengthBeam = 10;
beamType = 'GeometricallyExact';
numberOfElements = 8;
EI = 500;
EA = 10000;
nGP = 2;
rhoA = 0;
rhoI = 0;
kGA = EA;
allTimesteps = [5, 10, 50];

for index = 1:numel(allTimesteps)
    % setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = mfilename;
    setupObject.saveObject.saveData = false;
    setupObject.totalTimeSteps = allTimesteps(index);
    setupObject.totalTime = 1;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.newton.tolerance = 1e-11;
    setupObject.integrator = 'Midpoint';
    setupObject.plotObject.flag = true;
    setupObject.plotObject.postPlotType = "stress";
    setupObject.plotObject.stress.component = 3; %select 1 for bending moment and 2 for shear force
    setupObject.plotObject.view = [0,90];
    setupObject.plotObject.setXminXmax([-2;-0.1;-1],[10;8;1]);
    setupObject.plotObject.border = 0.0;
    setupObject.plotObject.plotInitialConfig = true;
    setupObject.plotObject.showGrid = true;
    setupObject.plotObject.colorscheme = 'viridis';
    
    dofObject = dofClass; % required object for dof and object handling
    dofObject.postDataObject.storeStateFlag = true;
    
    %% continuum Objects
    beamObject = beamClass(dofObject,2);
    beamObject.materialObject.name = "Hooke";
    beamObject.theory = beamType;
    beamObject.elementDisplacementType = 'mixedPH';
    beamObject.materialObject.EI = EI;
    beamObject.materialObject.EA = EA;
    beamObject.materialObject.kGA = kGA;
    beamObject.materialObject.rhoI = rhoI;
    beamObject.materialObject.rhoA = rhoA;
    beamObject.materialObject.rho = 1;
    
    [nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
    nodes = nodes + 1 / 2 * lengthBeam;
    beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];
    
    beamObject.dimension = 1;
    beamObject.shapeFunctionObject.order = 1;
    beamObject.shapeFunctionObject.numberOfGausspoints = nGP;
    
    % dirichlet boundaries
    dirichletObject = dirichletClass(dofObject);
    dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
    dirichletObject.nodalDof = [1, 2, 3];
    dirichletObject.masterObject = beamObject;
    
    % neumann boundaries
    Q = [0;0;2*EI*pi/lengthBeam];
    nodalLoadObject = nodalLoadClass(dofObject);
    nodalLoadObject.masterObject = beamObject;
    nodalLoadObject.loadVector = Q;
    nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);
    nodalLoadObject.timeFunction = str2func('@(t) t');
    
    %% solver
    beamObject.numericalTangentObject.computeNumericalTangent = true;
    beamObject.numericalTangentObject.showDifferences = false;
    dofObject = runNewton(setupObject, dofObject);
    
    %matlab2tikz('height', '\figH', 'width', '\figW', 'filename', ['snapshots_rollup_', num2str(allTimesteps(index)), '.tikz'], 'showInfo', false, 'floatformat', '%.4g')
    
end