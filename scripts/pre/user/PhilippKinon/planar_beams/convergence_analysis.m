clearvars;
%% Aufgabe 13 a) - Grundlagen Finite Elemente

lengthBeam = 100;
widthToLength = 1 / 50;

beamType = 'Timoshenko';
dispType = {'displacement', 'displacement', 'mixedPH', 'mixedPH'};
nGP = [2,1,2,1];
numberOfElements = [1,2,4,8,16,32];
displacementLastNode = zeros(numel(dispType),numel(numberOfElements));

for formulation = 1:numel(dispType)
    for nel = 1:numel(numberOfElements)
        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = mfilename;
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.newton.tolerance = 1e-7;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass; % required object for dof and object handling
        
        %% continuum Objects
        beamObject = beamClass(dofObject);
        beamObject.theory = beamType;
        beamObject.materialObject.name = 'Hooke';
        beamObject.elementDisplacementType = dispType{formulation}; 
        widthBeam = widthToLength * lengthBeam;
        beamObject.materialObject.E = 2.1 * 10^6;
        nu = 0.3;
        beamObject.materialObject.G = beamObject.materialObject.E/(2*(1+nu));
        beamObject.materialObject.I = widthBeam^4 / 12;
        beamObject.materialObject.A = widthBeam^2;
        beamObject.materialObject.rho = 0;
        beamObject.materialObject.shearCorrectionCoefficient = 5/6;
        
        [nodes, beamObject.meshObject.edof, edofNeumann] = meshOneDimensional(lengthBeam, numberOfElements(nel), 1);
        nodes = nodes + 1 / 2 * lengthBeam;
        beamObject.meshObject.nodes = [nodes, zeros(size(nodes))];
        
        beamObject.dimension = 1;
        beamObject.shapeFunctionObject.order = 1;
        beamObject.shapeFunctionObject.numberOfGausspoints = nGP(formulation);
        
        % dirichlet boundaries
        dirichletObject = dirichletClass(dofObject);
        dirichletObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == 0);
        dirichletObject.nodalDof = [1, 2];
        dirichletObject.masterObject = beamObject;
        
        % neumann boundaries
        Q = [100;0];
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = beamObject;
        nodalLoadObject.loadVector = Q;
        nodalLoadObject.nodeList = find(beamObject.meshObject.nodes(:, 1) == lengthBeam);
        
        %% solver
        dofObject = runNewton(setupObject, dofObject);
        displacementLastNode(formulation,nel) = dofObject.listContinuumObjects{1}.qN1(end, 1);
    end
end

% Compute increment
figure()
plot(numberOfElements, displacementLastNode(1,:), numberOfElements, displacementLastNode(2,:), numberOfElements, displacementLastNode(3,:))
lgd = legend('dispnGP2', 'dispnGP1', 'mixedPH');
lgd.Location = "southeast";
xlabel('nel')
ylabel('wendpoint')
grid on
matlab2tikz('height', '\figH', 'width', '\figW', 'filename', [mfilename,num2str(widthToLength),'.tikz'], 'showInfo', false, 'floatformat', '%.6g')