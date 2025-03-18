% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'easPetrovGalerkin'}; % {'disp', 'easPetrovGalerkin', 'dispPetrovGalerkinZhao'};
%% Script
resultArray = zeros(size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)

        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp('=======================================');

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'temp';
        setupObject.plotObject.stress.component = -1;
        setupObject.plotObject.view = 2;
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-3;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
        lengthX = 1;
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
                solidThermoObject.materialObject.name = 'HookeESZ';
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
        
        meshType = 'distorted4';
        if strcmp(meshType, 'regular')
            [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 20, 20, order, serendipity);
        elseif strcmp(meshType, 'distorted4')
            [solidThermoObject.meshObject.nodes,solidThermoObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 2, 2, order, serendipity);
            nodeToDistort = find(solidThermoObject.meshObject.nodes(:,1) == 0 & solidThermoObject.meshObject.nodes(:,2) == 0);
            solidThermoObject.meshObject.nodes(nodeToDistort, :) = solidThermoObject.meshObject.nodes(nodeToDistort, :) + lengthX/10;
        end
        solidThermoObject.meshObject.nodes = adaptMeshToCornerNodes(solidThermoObject.meshObject.nodes, solidThermoObject.meshObject.edof, 2);
        solidThermoObject.meshObject.nodes(:,1) = solidThermoObject.meshObject.nodes(:,1) + lengthX/2;
        solidThermoObject.meshObject.nodes(:,2) = solidThermoObject.meshObject.nodes(:,2) + lengthY/2;

        solidThermoObject.meshObject.nodes = [solidThermoObject.meshObject.nodes, 293.15*ones(size(solidThermoObject.meshObject.nodes,1), 1)];

        solidThermoObject.materialObject.rho = 0;
        solidThermoObject.materialObject.E = 1e7;
        solidThermoObject.materialObject.nu = 0.0;
        solidThermoObject.materialObject.alphaT = 1e-5;
        solidThermoObject.materialObject.k = 0.1;
        solidThermoObject.dimension = 2;
        solidThermoObject.shapeFunctionObject.order = order;
        solidThermoObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidThermoObject.meshObject.nodes(:,1) == 0);
        dirichletBoundary.nodalDof = 1:2;
        dirichletBoundary.masterObject = solidThermoObject;

        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidThermoObject.meshObject.nodes(:,1) == 0);
        dirichletBoundary2.nodalDof = 3;
        dirichletBoundary2.masterObject = solidThermoObject;
        dirichletBoundary2.timeFunction = @(t, XYZ) 293.15;
        
        % Neumann boundary
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.loadPhysics = 'thermal';
        neumannObject.masterObject = solidThermoObject;
        neumannObject.loadVector = 1; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX2;
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        nodeToMeasure = find(solidThermoObject.meshObject.nodes(:,1) == lengthX & solidThermoObject.meshObject.nodes(:,2) == 0);
        endDisplacement = solidThermoObject.qN1(nodeToMeasure, 2) - solidThermoObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        disp(['End temperature: ', num2str(solidThermoObject.qN1(nodeToMeasure, 3))])
        resultArray(ii) = endDisplacement;
        
end