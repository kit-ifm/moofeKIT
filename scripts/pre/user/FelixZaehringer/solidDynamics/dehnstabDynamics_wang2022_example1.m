%% Parameters
lengthX = 5;
lengthY = 0.5;
s=0.15; % 0.15
% s=0

% class = 'solid'; % solid, solidVelocity
class = 'solidVelocity'; % solid, solidVelocity

timeStepVector=[4, 8, 16, 20, 32, 40, 48, 60, 64, 80, 96, 100, 128, 144, 160, 176, 192, 200];
timeStepSizeVector = zeros(size(timeStepVector));
errorVector = zeros(size(timeStepVector));
for jj=1:length(timeStepVector)
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'dehnstabDynamics_wang2022_example1';
    setupObject.saveObject.saveData = false;
    setupObject.totalTime = 0.2;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.view = 2;
    setupObject.plotObject.postPlotType = 'zero';
    setupObject.plotObject.stress.component = 11;
    setupObject.newton.tolerance = 1e-2;
    % setupObject.integrator = 'Endpoint';
    setupObject.integrator = 'Midpoint';
    % setupObject.integrator = 'ExplicitEuler';
    % setupObject.integrator = 'Newmark';
    % setupObject.newmarkGamma         = 0.5;
    % setupObject.newmarkBeta          = 0.25;
    switch setupObject.integrator
        case {'Endpoint', 'Newmark', 'Midpoint'}
            setupObject.totalTimeSteps = timeStepVector(jj); %2000
        case {'ExplicitEuler'}
            setupObject.totalTimeSteps = 1000;
            %courant: dt<=4.5e-4 (12000 steps)
        otherwise
            error('not implemented')
    end

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    if strcmpi(class, 'solid')
        solidObject = solidClass(dofObject);
    elseif strcmpi(class, 'solidVelocity')
        solidObject = solidVelocityClass(dofObject);
        setupObject.disableTimeIntegration = true;
    end
    solidObject.dimension = 2;
    % solidObject.elementDisplacementType = 'mixed';
    % solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
    % solidObject.elementNameAdditionalSpecification = 'NewConstraint';
    solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinNewConstraint';
    solidObject.numericalTangentObject.computeNumericalTangent = false;
    solidObject.numericalTangentObject.showDifferences = false;
    % solidObject.mixedFEObject.typeShapeFunction = 1;
    % solidObject.mixedFEObject.typeShapeFunctionData = 16;
    solidObject.mixedFEObject.condensation = false;
    order = 2;
    serendipity = true;
    solidObject.shapeFunctionObject.order = order;
    solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
    solidObject.materialObject.name = 'HookeESZ';
    if strcmpi(class, 'solid')
        solidObject.materialObject.rho = 100;
    elseif strcmpi(class, 'solidVelocity')
        solidObject.materialObject.rho = 100;
        solidObject.materialObject.rhoNew = 100;
    end

    solidObject.materialObject.E = 1e6;
    solidObject.materialObject.nu = 0;
    v0 = 10;

    meshType = 'distorted';

    if strcmpi(meshType, 'simple')
        [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 20, 1, order, serendipity);
        solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
    elseif strcmpi(meshType, 'distorted')
        [nodes, solidObject.meshObject.edof, edofBoundary] = meshRectangle(lengthX, lengthY, 10, 1, order, serendipity);
        nodesToDistortTop = find(nodes(:, 1) ~= -lengthX/2 & nodes(:, 1) ~= lengthX/2 & nodes(:,2) == lengthY/2);
        nodesToDistortBottom = find(nodes(:, 1) ~= -lengthX/2 & nodes(:, 1) ~= lengthX/2 & nodes(:,2) == -lengthY/2);
        nodes(nodesToDistortTop(2:4:end), 1) = nodes(nodesToDistortTop(2:4:end), 1) + s;
        nodes(nodesToDistortTop(4:4:end), 1) = nodes(nodesToDistortTop(4:4:end), 1) - s;
        nodes(nodesToDistortBottom(2:4:end), 1) = nodes(nodesToDistortBottom(2:4:end), 1) - s;
        nodes(nodesToDistortBottom(4:4:end), 1) = nodes(nodesToDistortBottom(4:4:end), 1) + s;
        nodes = adaptMeshToCornerNodes(nodes, solidObject.meshObject.edof, 2);
        solidObject.meshObject.nodes(:,1) = nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = nodes(:,2) + lengthY/2;
    end

    if strcmpi(class, 'solidVelocity')
        solidObject.meshObject.nodes = [solidObject.meshObject.nodes, [v0*ones(size(solidObject.meshObject.nodes, 1), 1), zeros(size(solidObject.meshObject.nodes, 1), 1)]];
    end

    % dirchlet boundary conditions
    dirichletObject1 = dirichletClass(dofObject);
    dirichletObject1.masterObject = solidObject;
    dirichletObject1.nodeList = find(solidObject.meshObject.nodes(:, 1) == lengthX);
    dirichletObject1.timeFunction = @(t, XYZ) XYZ;
    dirichletObject1.nodalDof = 1;

    dirichletObject2 = dirichletClass(dofObject);
    dirichletObject2.masterObject = solidObject;
    dirichletObject2.nodeList = find(solidObject.meshObject.nodes(:, 1) == lengthX);
    dirichletObject2.timeFunction = @(t, XYZ) XYZ;
    dirichletObject2.nodalDof = 2;

    if strcmpi(class, 'solidVelocity')
        dirichletObject3 = dirichletClass(dofObject);
        dirichletObject3.masterObject = solidObject;
        dirichletObject3.nodeList = find(solidObject.meshObject.nodes(:, 1) == lengthX);
        dirichletObject3.nodalDof = 3;

        dirichletObject4 = dirichletClass(dofObject);
        dirichletObject4.masterObject = solidObject;
        dirichletObject4.nodeList = find(solidObject.meshObject.nodes(:, 1) == lengthX);
        dirichletObject4.nodalDof = 4;

        % set initial velocity at bearings to zero => necessary to show energy conservation
        solidObject.meshObject.nodes(dirichletObject3.nodeList, 3) = 0;
    end

    %% solver
    dofObject = runNewton(setupObject, dofObject);

    % FEM solution
    nodeToTrack = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
    xDisplacementNodeToTrack = solidObject.qN1(nodeToTrack, 1) - solidObject.qR(nodeToTrack, 1);
    timeStepSizeVector(jj) = setupObject.timeStepSize;
    disp(['Time step size: ', num2str(setupObject.timeStepSize)]);

    dispError = getElementData(dofObject.postDataObject, dofObject, setupObject, 'dispError');
    errorVector(jj) = dispError(end);

    % disp(['Error: ', num2str(abs(xDisplacementNodeToTrack - (-0.199999999954504))/abs(-0.199999999954504), '%4.16f')]);

end

figure;
log10log10(timeStepSizeVector, errorVector);

% Lineare Regression im log-log-Raum
p = polyfit(log(timeStepSizeVector(1:4)), log(errorVector(1:4)), 1);

% Die Steigung der Ausgleichsgeraden
steigung = p(1)


% exact solution
% xAna = 0;
% xEval = 0; % lengthX/2
% tEval=0:0.0001:setupObject.totalTime;
% tEval=0.32;
% for n=1:1000
%     c = sqrt(solidObject.materialObject.E/solidObject.materialObject.rho);
%     omega_n = (2*n-1)*pi*c/(2*lengthX);
%     A_n = 8*lengthX*v0*(-1)^(n+1)/(((2*n-1)*pi)^2*c);
%     xAna = xAna + A_n*sin(omega_n*tEval)*cos((2*n-1)*pi*xEval/(2*lengthX));
% end
%
% figure;
% plot(tEval(1:20:end), xAna(1:20:end));
% xlim([0, setupObject.totalTime]);
% ylim([-0.8, 0.8]);

%% postprocessing - energy
% timeVector = getTime(dofObject.postDataObject, setupObject);
% if strcmpi(class, 'solid')
%     kineticEnergy = getKineticEnergy(dofObject.postDataObject, setupObject);
% elseif strcmpi(class, 'solidVelocity')
%     kineticEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'kineticEnergy');
% end
% strainEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'strainEnergy');
% externalEnergy = getElementData(dofObject.postDataObject, dofObject, setupObject, 'externalEnergy');
% totalEnergy = kineticEnergy + strainEnergy;
% totalEnergyDifference = totalEnergy(2:end) - totalEnergy(1:end-1);
% figure;
% plot(timeVector(1:2:end), kineticEnergy(1:2:end)+strainEnergy(1:2:end));
% figure;
% plot(timeVector(1:2:end), strainEnergy(1:2:end));
% figure;
% plot(timeVector(1:2:end-1), totalEnergyDifference(1:2:end));

%% postprocessing
% [setupObjectCell, dofObjectCell] = loadObjectsFromMatFile('dehnstabDynamics_wang2022_example1', 'all');
% nodeToTrack = find(solidObject.meshObject.nodes(:, 1) == 0 & solidObject.meshObject.nodes(:, 2) == lengthY/2);
% xDisplacementNodeToTrack = zeros(size(dofObjectCell, 1), 1);
% timeVector = zeros(size(dofObjectCell, 1), 1);
% for i=1:size(dofObjectCell, 1)
%     setupObjectOfTimestep = setupObjectCell{i};
%     dofObjectOfTimestep = dofObjectCell{i};
%     solidObjectOfTimestep = dofObjectOfTimestep;
%     xDisplacementNodeToTrack(i) = solidObjectOfTimestep(nodeToTrack, 1) - solidObject.qR(nodeToTrack, 1);
%     timeVector(i) = setupObjectOfTimestep;
% end
% figure;
% plot([0; timeVector(1:2:end)], [0; xDisplacementNodeToTrack(1:2:end)]);
% xlim([0, setupObject.totalTime]);
% ylim([-0.8, 0.8]);
