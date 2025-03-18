%% Start parameters
elementsToTest = {'displacementQ8', 'displacementPetrovGalerkinQ8'};
numberStepsMeshDistortion = 9;
maxDistortion = 4.5;
testBothDirections = false;

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
        setupObject.plotObject.stress.component = 22;
        setupObject.plotObject.view = 2;
        setupObject.integrator = 'Endpoint';
        setupObject.newton.tolerance = 1e-4;

        dofObject = dofClass; % required object for dof and object handling

        %% continuum Objects
        meshType = 'B'; % A or B (regular/distorted)
        lengthX = 10;
        lengthY = 1; % 0.1
        axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
        axisymmetricSolidObject.dimension = 2;

        % axisymmetricSolidObject.numericalTangentObject.computeNumericalTangent = true;
        % axisymmetricSolidObject.numericalTangentObject.showDifferences = true;
        % axisymmetricSolidObject.numericalTangentObject.type = 'complex';

        numberOfElementsX = 2;
        numberOfElementsY = 1;

        switch elementsToTest{ii}
            case 'displacementQ8'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'displacementQ4'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 1;
                serendipity = false;
                numberOfElementsX = 80;
                numberOfElementsY = 10;
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
                numberOfElementsX = 4;
                numberOfElementsY = 2;
            case 'easKasperTaylor'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'eas';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'KasperTaylor';
                axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 2;
                axisymmetricSolidObject.mixedFEObject.condensation = true;
                order = 1;
                serendipity = false;
                numberOfElementsX = 4;
                numberOfElementsY = 2;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        innerRadius = 0;
        axisymmetricSolidObject.shapeFunctionObject.order = order;
        axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

        if strcmp(meshType, 'A')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
        elseif strcmp(meshType, 'B')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
            nodesToDistortBottom = find(axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
            nodesToDistortTop = find(axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius+lengthX & axisymmetricSolidObject.meshObject.nodes(:, 1) ~= innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
            nodeMiddle = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            if order == 2
                nodesToDistortBottom = nodesToDistortBottom(2:2:end);
                nodesToDistortTop = nodesToDistortTop(2:2:end);
            end
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) + s;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) - s;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) - s;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) + s;
            axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) = axisymmetricSolidObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) + s;
        elseif strcmp(meshType, 'BCurvedEdge')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX / 2 + innerRadius;
            nodeToDistort = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius+lengthX/2 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) + s;
        elseif strcmp(meshType, 'C')
            [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity);
            axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + innerRadius;
            axisymmetricSolidObject.meshObject.nodes(:, 2) = axisymmetricSolidObject.meshObject.nodes(:, 2) - lengthY / 2;
        end
        axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);
%         axisymmetricSolidObject.meshObject.nodes(nodeMiddle, 1) = axisymmetricSolidObject.meshObject.nodes(nodeMiddle, 1) + s/2;

        axisymmetricSolidObject.materialObject.rho = 0;
        axisymmetricSolidObject.materialObject.E = 1e7;
        axisymmetricSolidObject.materialObject.nu = 0.3;

        % Dirichlet boundary
        boundaryCondition = 'doubleSupport'; % clamped or simplySupported or doubleSupport
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
            dirichletBoundary2 = dirichletClass(dofObject);
            dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
            dirichletBoundary2.nodalDof = 2;
            dirichletBoundary2.timeFunction = @(t, XYZ) -lengthY / 2;
            dirichletBoundary2.masterObject = axisymmetricSolidObject;
        elseif strcmp(boundaryCondition, 'doubleSupport')
            dirichletBoundary = dirichletClass(dofObject);
            dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            dirichletBoundary.nodalDof = 2;
            dirichletBoundary.timeFunction = @(t, XYZ) 0;
            dirichletBoundary.masterObject = axisymmetricSolidObject;

            dirichletBoundary2 = dirichletClass(dofObject);
            dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == lengthX & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
            dirichletBoundary2.nodalDof = 1;
            dirichletBoundary2.timeFunction = @(t, XYZ) lengthX;
            dirichletBoundary2.masterObject = axisymmetricSolidObject;
            %
            %     dirichletBoundary2 = dirichletClass(dofObject);
            %     dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0);
            %     dirichletBoundary2.nodalDof = 1;
            %     dirichletBoundary2.masterObject = axisymmetricSolidObject;
        end

        % Neumann boundary
        % B matrix plate
        B = axisymmetricSolidObject.materialObject.E*lengthY^3/(12*(1-axisymmetricSolidObject.materialObject.nu^2));

        loadCondition = 'endMoment'; % lineLoad or endMoment
        if strcmp(loadCondition, 'lineLoad')
            neumannObject = neumannClass(dofObject);
            neumannObject.loadGeometry = 'line';
            neumannObject.masterObject = axisymmetricSolidObject;
            neumannObject.loadVector = [0; -10.24]; % Komponenten des Kraftvektors anpassen
            neumannObject.meshObject.edof = bounEdof.SY2;
        elseif strcmp(loadCondition, 'endMoment')
            M0 = 2*(1+axisymmetricSolidObject.materialObject.nu)*B/(lengthX^2);
            neumannObject = neumannClass(dofObject);
            neumannObject.loadGeometry = 'line';
            neumannObject.masterObject = axisymmetricSolidObject;
            neumannObject.loadVectorFunction = @(XYZ) [-12/(lengthY)^3*M0*XYZ(2); 0]; % Komponenten des Kraftvektors anpassen
            neumannObject.meshObject.edof = bounEdof.SX2;
        end

        %% solver
        dofObject = runNewton(setupObject, dofObject);
        nodeToMeasure = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == innerRadius & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
        endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 2) - axisymmetricSolidObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = -endDisplacement;
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
