% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'dispQ9'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
numberStepsMeshDistortion = 1;
maxDistortion = 25;

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
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 11;
        setupObject.plotObject.view = 2;
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
        lengthX = 100;
        lengthY = 10;
        solidObject = solidClass(dofObject);

        switch elementsToTest{ii}
            case 'disp'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeESZ';
            case 'dispQ4'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
            case 'dispQ9'
                order = 2;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
            case 'eas'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'eas';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkin'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkin';
            case 'dispPetrovGalerkinXie'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkinCen'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'displacement';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinCenEtAl2015';
            case 'easPetrovGalerkin'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                solidObject.elementDisplacementType = 'eas';
                solidObject.mixedFEObject.condensation = false;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        [solidObject.meshObject.nodes,solidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, 2, 1, order, serendipity);
        nodeToDistort1 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == -lengthY/2);
        nodeToDistort2 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY/2);
        solidObject.meshObject.nodes(nodeToDistort1, 1) = solidObject.meshObject.nodes(nodeToDistort1, 1) - s;
        solidObject.meshObject.nodes(nodeToDistort2, 1) = solidObject.meshObject.nodes(nodeToDistort2, 1) + s;
        solidObject.meshObject.nodes = adaptMeshToCornerNodes(solidObject.meshObject.nodes, solidObject.meshObject.edof, 2);
        solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2) + lengthY/2;

        %plateObject.materialObject.name = 'Hooke';
        solidObject.materialObject.rho = 0;
        solidObject.materialObject.E = 1e7;
        solidObject.materialObject.nu = 0.3;
        solidObject.dimension = 2;
        solidObject.shapeFunctionObject.order = order;
        solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary.nodalDof = 1:2;
        dirichletBoundary.masterObject = solidObject;

        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY);
        dirichletBoundary2.nodalDof = 1;
        dirichletBoundary2.masterObject = solidObject;
        

        % Neumann boundary
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) [0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX2;
    
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) -[0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX1;
    
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) [240*XYZ(2)/lengthY-120; 0]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX1;
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        nodeToMeasure = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == 0);
        endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement;
        
    end
end