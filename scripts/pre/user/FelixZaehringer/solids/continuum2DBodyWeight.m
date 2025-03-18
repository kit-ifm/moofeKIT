% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'disp', 'dispPetrovGalerkin'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
numberStepsMeshDistortion = 1;
maxDistortion = 9;

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
        setupObject.plotObject.stress.component = -1;
        setupObject.plotObject.view = 2;
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-2;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
        lengthX = 10;
        lengthY = 20;
        solidObject = solidClass(dofObject);

        switch elementsToTest{ii}
            case 'disp'
                order = 2;
                serendipity = true;
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
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkinXie'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easPetrovGalerkin'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'easPetrovGalerkin';
                solidObject.mixedFEObject.condensation = false;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshRectangle(lengthX, lengthY, 10, 20, order, serendipity);
        nodeToDistort1 = find(solidObject.meshObject.nodes(:,1) == -lengthX/2 & solidObject.meshObject.nodes(:,2) == 0);
        nodeToDistort2 = find(solidObject.meshObject.nodes(:,1) == lengthX/2 & solidObject.meshObject.nodes(:,2) == 0);
        solidObject.meshObject.nodes(nodeToDistort1, 2) = solidObject.meshObject.nodes(nodeToDistort1, 2) - s;
        solidObject.meshObject.nodes(nodeToDistort2, 2) = solidObject.meshObject.nodes(nodeToDistort2, 2) + s;
        solidObject.meshObject.nodes = adaptMeshToCornerNodes(solidObject.meshObject.nodes, solidObject.meshObject.edof, 2);
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
        dirichletBoundary.nodalDof = 1;
        dirichletBoundary.masterObject = solidObject;

        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary2.nodalDof = 2;
        dirichletBoundary2.masterObject = solidObject;
        
        % Body Forces
        px = 100*50;
        bodyForceObject = bodyForceClass(dofObject);
        bodyForceObject.masterObject = solidObject;
        bodyForceObject.dimension = solidObject.dimension;
        bodyForceObject.typeOfLoad = 'deadLoad';
        bodyForceObject.shapeFunctionObject.numberOfGausspoints = solidObject.shapeFunctionObject.numberOfGausspoints;
        bodyForceObject.meshObject.edof = solidObject.meshObject.edof;
        bodyForceObject.timeFunction = @(t) 1;
        bodyForceObject.loadFunction = [0; -px];


%         % Test with gleichlast
%         nodalForceObject = nodalForceClass(dofObject);
%         nodalForceObject.dimension = 2;                                % Dimension des zu betrachtenden Körpers angeben
%         nodalForceObject.masterObject = solidObject;
%         nodalForceObject.forceVector = [0; 200];                         % Komponenten des Kraftvektors anpassen
%         nodalForceObject.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == lengthY);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        nodeToMeasure = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY);
        endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement;
        
    end
end