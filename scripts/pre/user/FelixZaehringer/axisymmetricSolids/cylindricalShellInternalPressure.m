%% Start parameters
elementsToTest = {'displacementQ8', 'displacementPetrovGalerkinQ8'}; % , 'displacementQ8', 'displacementPetrovGalerkinQ8'
numberStepsMeshDistortion = 20;
maxDistortion = 4.8;
testBothDirections = true;

%% Script
totalNumberOfSteps = numberStepsMeshDistortion + 1;
if testBothDirections
    totalNumberOfSteps = 2 * numberStepsMeshDistortion + 1;
end
resultArray = zeros(totalNumberOfSteps, size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)
    for jj = 1:totalNumberOfSteps
        s = (jj - 1) * maxDistortion / numberStepsMeshDistortion;
        if testBothDirections
            s = (jj - 1 - numberStepsMeshDistortion) * maxDistortion / numberStepsMeshDistortion;
        end

        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp(['Current s: ', num2str(s)]);
        disp('=======================================');

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 11;
        setupObject.plotObject.view = 2;
        setupObject.integrator = 'Endpoint';
        setupObject.newton.tolerance = 1e-4;

        dofObject = dofClass; % required object for dof and object handling

        %% continuum Objects
        meshType = 'BCurvedEdge'; % A or B (regular/distorted)
        lengthX = 1;
        lengthY = 10; % 0.1
        axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
        axisymmetricSolidObject.dimension = 2;

        % axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
        % axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
        % axisymmetricSolidObject.numericalTangentObject.type = 'complex';

        switch elementsToTest{ii}
            case 'displacementQ8'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'displacementQ4'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinQ4'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'displacementQ9'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = false;
            case 'displacementPetrovGalerkinQ8'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'displacement';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'displacementPetrovGalerkinQ9'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'displacement';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = false;
            case 'easSimoRifai'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'eas';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'SimoRifai';
                axisymmetricSolidObject.mixedFEObject.condensation = true;
                axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 5;
                order = 1;
                serendipity = false;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        numberOfElementsX = 1;
        numberOfElementsY = 2;
        innerRadius = 4.5;
        axisymmetricSolidObject.shapeFunctionObject.order = order;
        axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

        if strcmp(meshType, 'A')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
            axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) + lengthY / 2;
        elseif strcmp(meshType, 'B')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
            nodesToDistortLeft = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) ~= lengthY/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) ~= -lengthY/2);
            nodesToDistortRight = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) ~= lengthY/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) ~= -lengthY/2);
            nodeMiddle = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX/2);
%             if order == 2
%                 nodesToDistortLeft = nodesToDistortLeft(2:2:end);
%                 nodesToDistortRight = nodesToDistortRight(2:2:end);
%             end
%             axisymmetricSolidObject.meshObject.nodes(nodesToDistortLeft(1:2:end), 2) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortLeft(1:2:end), 2) + s/2;
%             axisymmetricSolidObject.meshObject.nodes(nodesToDistortLeft(2:2:end), 2) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortLeft(2:2:end), 2) + s;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortRight(1:2:end), 2) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortRight(1:2:end), 2) - s/2;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortRight(2:2:end), 2) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortRight(2:2:end), 2) - s;
            axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) + lengthY / 2;
        elseif strcmp(meshType, 'BCurvedEdge')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
            nodeToDistort = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) + s;
            axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) + lengthY / 2;
        elseif strcmp(meshType, 'C')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + innerRadius;
            axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) - lengthY / 2;
        end
%         axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);
%         axisymmetricSolidObject.meshObject.nodes(nodeMiddle, 1) = axisymmetricSolidObject.meshObject.nodes(nodeMiddle, 1) + s/2;
        

        axisymmetricSolidObject.materialObject.rho = 0;
        axisymmetricSolidObject.materialObject.E = 10.92e5;
        axisymmetricSolidObject.materialObject.nu = 0.3;

        % Dirichlet boundary
        boundaryCondition = 'simplySupported'; % clamped or simplySupported or doubleSupport
        if strcmp(boundaryCondition, 'clamped')
            dirichletBoundary = dirichletClass(dofObject);
            dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX);
            dirichletBoundary.nodalDof = 1;
            dirichletBoundary.timeFunction = @(t, XYZ) lengthX;
            dirichletBoundary.masterObject = axisymmetricSolidObject;

            dirichletBoundary2 = dirichletClass(dofObject);
            dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
            dirichletBoundary2.nodalDof = 2;
            dirichletBoundary2.timeFunction = @(t, XYZ) -lengthY / 2;
            dirichletBoundary2.masterObject = axisymmetricSolidObject;

            dirichletBoundary3 = dirichletClass(dofObject);
            dirichletBoundary3.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
            dirichletBoundary3.nodalDof = 2;
            dirichletBoundary3.timeFunction = @(t, XYZ) lengthY / 2;
            dirichletBoundary3.masterObject = axisymmetricSolidObject;
        elseif strcmp(boundaryCondition, 'simplySupported')
            dirichletBoundary = dirichletClass(dofObject);
            dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius + lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            dirichletBoundary.nodalDof = 2;
            dirichletBoundary.timeFunction = @(t, XYZ) 0;
            dirichletBoundary.masterObject = axisymmetricSolidObject;
        elseif strcmp(boundaryCondition, 'doubleSupport')
            dirichletBoundary = dirichletClass(dofObject);
            dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius + lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
%             dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            dirichletBoundary.nodalDof = 2;
            dirichletBoundary.timeFunction = @(t, XYZ) 0;
            dirichletBoundary.masterObject = axisymmetricSolidObject;

            dirichletBoundary2 = dirichletClass(dofObject);
            dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius + lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
%             dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            dirichletBoundary2.nodalDof = 1;
            dirichletBoundary2.timeFunction = @(t, XYZ) XYZ;
            dirichletBoundary2.masterObject = axisymmetricSolidObject;
            %
            %     dirichletBoundary2 = dirichletClass(dofObject);
            %     dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0);
            %     dirichletBoundary2.nodalDof = 1;
            %     dirichletBoundary2.masterObject = axisymmetricSolidObject;
        end

        % Neumann boundary
        loadCondition = 'lineLoadFunction'; % lineLoad or endMoment
        if strcmp(loadCondition, 'lineLoad')
            neumannObject = neumannClass(dofObject);
            neumannObject.loadGeometry = 'line';
            neumannObject.masterObject = axisymmetricSolidObject;
            neumannObject.loadVector = [1e5; 0]; % Komponenten des Kraftvektors anpassen
            neumannObject.meshObject.edof = bounEdof.SX1;
        elseif strcmp(loadCondition, 'lineLoadFunction')
            gamma = 5000;
            neumannObject = neumannClass(dofObject);
            neumannObject.loadGeometry = 'line';
            neumannObject.masterObject = axisymmetricSolidObject;
            neumannObject.loadVectorFunction = @(XYZ) [gamma*lengthY-gamma*XYZ(2); 0]; % Komponenten des Kraftvektors anpassen
            neumannObject.meshObject.edof = bounEdof.SX1;
        elseif strcmp(loadCondition, 'endMoment')
            M0 = 20800/8*1e3;
            nodalLoadObject = nodalLoadClass(dofObject);
            nodalLoadObject.masterObject = axisymmetricSolidObject;
            nodalLoadObject.loadVector = M0 / lengthY * [0; 1] * 2 * pi * lengthX;
            nodalLoadObject.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY);

            nodalLoadObject2 = nodalLoadClass(dofObject);
            nodalLoadObject2.masterObject = axisymmetricSolidObject;
            nodalLoadObject2.loadVector = M0 / lengthY * [0; -1] * 2 * pi * lengthX;
            nodalLoadObject2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY);
        end

        % Analytical solution
        B = axisymmetricSolidObject.materialObject.E*lengthY^3/(12*(1-axisymmetricSolidObject.materialObject.nu^2));
%         w = (innerRadius+lengthX/2)/(axisymmetricSolidObject.materialObject.E*lengthX)*(innerRadius+lengthX/2)*neumannObject.loadVector(1);
%         u = axisymmetricSolidObject.materialObject.nu*gamma*(innerRadius+lengthX/2)/(2*axisymmetricSolidObject.materialObject.E*lengthX)*(lengthY^2);
%         w = gamma*(innerRadius+lengthX/2)^2/(axisymmetricSolidObject.materialObject.E*lengthX)*lengthY;

        %% solver
        dofObject = runNewton(setupObject, dofObject);
        nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius + lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
        endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 1) - axisymmetricSolidObject.qN(nodeToMeasure, 1);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement;
    end
end
figure;
if testBothDirections
    xVals = -maxDistortion:maxDistortion / numberStepsMeshDistortion:maxDistortion;
else
    xVals = 0:maxDistortion / numberStepsMeshDistortion:maxDistortion;
end

for ii = 1:size(elementsToTest, 2)
    plot(xVals, resultArray(:, ii)');
    hold on;
end
xlim([min(xVals), max(xVals)]);
% ylim([min(min(resultArray)), max(max(resultArray))]);
% ylim([0.0, 0.9]);
xlabel('s');
ylabel('w_{max}');
legend(elementsToTest);

