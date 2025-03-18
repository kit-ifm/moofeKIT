% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'eas', 'easPetrovGalerkin'}; % {'dispQ4', 'easPetrovGalerkin', 'dispPetrovGalerkinZhao'};
numberStepsMeshDistortion = 100;
maxDistortion = 4.8;

%% Script
resultArray = zeros(numberStepsMeshDistortion, size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)
    for jj = 1:numberStepsMeshDistortion+1
        s = (jj - 1)*maxDistortion/numberStepsMeshDistortion;

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
        setupObject.plotObject.postPlotType = 'temp';
        setupObject.plotObject.stress.component = -1;
        setupObject.plotObject.view = 2;
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
        lengthX = 10;
        lengthY = 1;
        solidThermoObject = solidThermoClass(dofObject);

        switch elementsToTest{ii}
            case 'disp'
                order = 2;
                serendipity = true;
                solidThermoObject.materialObject.name = 'HookeESZ';
            case 'dispQ4'
                order = 1;
                serendipity = false;
                solidThermoObject.materialObject.name = 'HookeEVZ';
            case 'dispQ9'
                order = 2;
                serendipity = false;
                solidThermoObject.materialObject.name = 'HookeESZ';
            case 'eas'
                order = 1;
                serendipity = false;
                solidThermoObject.materialObject.name = 'HookeESZ';
                solidThermoObject.elementDisplacementType = 'eas';
                solidThermoObject.mixedFEObject.condensation = true;
                solidThermoObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkin'
                order = 2;
                serendipity = true;
                solidThermoObject.materialObject.name = 'HookeESZ';
                solidThermoObject.elementDisplacementType = 'displacementPetrovGalerkin';
            case 'dispPetrovGalerkinXie'
                order = 2;
                serendipity = true;
                solidThermoObject.materialObject.name = 'HookeESZ';
                solidThermoObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
                solidThermoObject.mixedFEObject.condensation = true;
                solidThermoObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkinZhao'
                order = 1;
                serendipity = false;
                solidThermoObject.materialObject.name = 'HookeESZ';
                solidThermoObject.elementDisplacementType = 'displacement';
                solidThermoObject.elementNameAdditionalSpecification = 'PetrovGalerkinZhaoEtAl2024';
            case 'easPetrovGalerkin'
                order = 1;
                serendipity = false;
                solidThermoObject.materialObject.name = 'HookeESZ';
                solidThermoObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                solidThermoObject.elementDisplacementType = 'eas';
                solidThermoObject.mixedFEObject.condensation = false;
                solidThermoObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        meshType = 'Rajendran'; % Zhao or Rajendran

        if strcmp(meshType, 'Rajendran')
            [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 2, 1, order, serendipity);
            nodeToDistort1 = find(solidThermoObject.meshObject.nodes(:,1) == 0 & solidThermoObject.meshObject.nodes(:,2) == -lengthY/2);
            nodeToDistort2 = find(solidThermoObject.meshObject.nodes(:,1) == 0 & solidThermoObject.meshObject.nodes(:,2) == lengthY/2);
            solidThermoObject.meshObject.nodes(nodeToDistort1, 1) = solidThermoObject.meshObject.nodes(nodeToDistort1, 1) - s;
            solidThermoObject.meshObject.nodes(nodeToDistort2, 1) = solidThermoObject.meshObject.nodes(nodeToDistort2, 1) + s;
            solidThermoObject.meshObject.nodes = adaptMeshToCornerNodes(solidThermoObject.meshObject.nodes, solidThermoObject.meshObject.edof, 2);
        elseif strcmp(meshType, 'ZhaoII')
            [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 6, 2, order, serendipity);
            nodesToDistortTop = find(solidThermoObject.meshObject.nodes(:,2) == lengthY/2 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            nodesToDistortMiddle = find(solidThermoObject.meshObject.nodes(:,2) == 0 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            nodesToDistortBottom = find(solidThermoObject.meshObject.nodes(:,2) == -lengthY/2 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            solidThermoObject.meshObject.nodes(nodesToDistortTop, 1) = solidThermoObject.meshObject.nodes(nodesToDistortTop, 1) + s;
            solidThermoObject.meshObject.nodes(nodesToDistortMiddle, 1) = solidThermoObject.meshObject.nodes(nodesToDistortMiddle, 1) - s;
            solidThermoObject.meshObject.nodes(nodesToDistortBottom, 1) = solidThermoObject.meshObject.nodes(nodesToDistortBottom, 1) + s;
            solidThermoObject.meshObject.nodes = adaptMeshToCornerNodes(solidThermoObject.meshObject.nodes, solidThermoObject.meshObject.edof, 2);
        elseif strcmp(meshType, 'ZhaoIII')
            [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 6, 2, order, serendipity);
            nodesToDistortTop = find(solidThermoObject.meshObject.nodes(:,2) == lengthY/2 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            nodesToDistortMiddle = find(solidThermoObject.meshObject.nodes(:,2) == 0 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            nodesToDistortBottom = find(solidThermoObject.meshObject.nodes(:,2) == -lengthY/2 & solidThermoObject.meshObject.nodes(:,1) ~= -lengthX/2 & solidThermoObject.meshObject.nodes(:,1) ~= lengthX/2);
            solidThermoObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortTop(1:2:end), 1) + s;
            solidThermoObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortTop(2:2:end), 1) - s;
            solidThermoObject.meshObject.nodes(nodesToDistortMiddle(1:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortMiddle(1:2:end), 1) - s;
            solidThermoObject.meshObject.nodes(nodesToDistortMiddle(2:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortMiddle(2:2:end), 1) + s;
            solidThermoObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortBottom(1:2:end), 1) + s;
            solidThermoObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) = solidThermoObject.meshObject.nodes(nodesToDistortBottom(2:2:end), 1) - s;
            solidThermoObject.meshObject.nodes = adaptMeshToCornerNodes(solidThermoObject.meshObject.nodes, solidThermoObject.meshObject.edof, 2);
        end
        
        solidThermoObject.meshObject.nodes(:,1) = solidThermoObject.meshObject.nodes(:,1) + lengthX/2;
        solidThermoObject.meshObject.nodes(:,2) = solidThermoObject.meshObject.nodes(:,2) + lengthY/2;

        solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, 293.15*ones(size(solidThermoObject.meshObject.nodes,1), 1)];

        %plateObject.materialObject.name = 'Hooke';
        solidThermoObject.materialObject.rho = 0;
        solidThermoObject.materialObject.E = 1e6;
        solidThermoObject.materialObject.nu = 0.0;
        solidThermoObject.materialObject.alphaT = 1e-5;
        solidThermoObject.materialObject.k = 0.1;
        solidThermoObject.dimension = 2;
        solidThermoObject.shapeFunctionObject.order = order;
        solidThermoObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidThermoObject.meshObject.nodes(:,1) == 0);
        dirichletBoundary.nodalDof = 1;
        dirichletBoundary.masterObject = solidThermoObject;

        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidThermoObject.meshObject.nodes(:,1) == 0 & solidThermoObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary2.timeFunction = @(t, XYZ) XYZ;
        dirichletBoundary2.nodalDof = 2;
        dirichletBoundary2.masterObject = solidThermoObject;

        % thermal Dirichlet boundary
        dirichletBoundary3 = dirichletClass(dofObject);
        dirichletBoundary3.nodeList = find(solidThermoObject.meshObject.nodes(:,2) == lengthY);
        dirichletBoundary3.nodalDof = 3;
        dirichletBoundary3.timeFunction = @(t, XYZ) 293.15 - 10;
        dirichletBoundary3.masterObject = solidThermoObject;

        dirichletBoundary4 = dirichletClass(dofObject);
        dirichletBoundary4.nodeList = find(solidThermoObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary4.nodalDof = 3;
        dirichletBoundary4.timeFunction = @(t, XYZ) 293.15 + 10;
        dirichletBoundary4.masterObject = solidThermoObject;
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        nodeToMeasure = find(solidThermoObject.meshObject.nodes(:,1) == lengthX & solidThermoObject.meshObject.nodes(:,2) == lengthY);
        endDisplacement = solidThermoObject.qN1(nodeToMeasure, 2) - solidThermoObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement*100;
        
    end
end

figure;
xVals = 0:maxDistortion / numberStepsMeshDistortion:maxDistortion;

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