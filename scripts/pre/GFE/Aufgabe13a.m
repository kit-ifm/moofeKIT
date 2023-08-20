%% Aufgabe 13 a) - Grundlagen Finite Elemente

lengthBeam = 100;
widthToLength = [1 / 50, 1 / 3];
beamTypes = {'Bernoulli', 'Bernoulli'};
displacementLastNodeArray = zeros(length(widthToLength), length(beamTypes));
numberOfElementsArray = [1, 5; 1, 50];

for i = 1:length(widthToLength)
    for j = 1:length(beamTypes)

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'GFEAufgabe13a';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';

        dofObject = dofClass; % required object for dof and object handling

        %% continuum Objects
        beamObject = beamClass(dofObject);
        beamObject.materialObject.name = beamTypes{j};
        beamObject.elementDisplacementType = 'beam';
        widthBeam = widthToLength(i) * lengthBeam;
        beamObject.materialObject.E = 2.1 * 10^6;
        nu = 0.3;
        beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
        beamObject.materialObject.I = widthBeam^4 / 12;
        beamObject.materialObject.A = widthBeam^2;
        beamObject.materialObject.rho = 0;


        numberOfElements = numberOfElementsArray(i, j);

        [nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements, 1);
        nodes = nodes + 1 / 2 * lengthBeam;
        beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];

        beamObject.dimension = 1;
        beamObject.shapeFunctionObject.order = 1;
        beamObject.shapeFunctionObject.numberOfGausspoints = 2;

        % dirichlet boundaries
        dirichletObject = dirichletClass(dofObject);
        dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
        dirichletObject.nodalDof = [1, 2];
        dirichletObject.masterObject = beamObject;

        % neumann boundaries
        Q = 100;
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = beamObject;
        nodalLoadObject.loadVector = -Q;
        nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);

        %% solver
        dofObject = runNewton(setupObject, dofObject);
        displacementLastNodeArray(i, j) = dofObject.listContinuumObjects{1}.qN1(end, 1);
    end
end