%% Start parameters
elementsToTest = {'dispPetrovGalerkinXie'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD'};
numberStepsCoordinateRotation = 3;
maxCoordinateRotation = 25;

%% Script
resultArray = zeros(numberStepsCoordinateRotation, size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)
    for jj = 1:numberStepsCoordinateRotation+1
        phi = (jj - 1)*maxCoordinateRotation/numberStepsCoordinateRotation;

        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp(['Current phi: ', num2str(s)]);
        disp('=======================================');

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-5;
        
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
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easPetrovGalerkin'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'easPetrovGalerkin';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

%         [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshRectangle(lengthX, lengthY, 2, 1, order, true);
%         nodeToDistort1 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == -lengthY/2);
%         nodeToDistort2 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY/2);
%         solidObject.meshObject.nodes(nodeToDistort1, 1) = solidObject.meshObject.nodes(nodeToDistort1, 1) - s;
%         solidObject.meshObject.nodes(nodeToDistort1-1, 1) = solidObject.meshObject.nodes(nodeToDistort1-1, 1) - s/2;
%         solidObject.meshObject.nodes(nodeToDistort1+1, 1) = solidObject.meshObject.nodes(nodeToDistort1+1, 1) - s/2;
%         solidObject.meshObject.nodes(nodeToDistort2, 1) = solidObject.meshObject.nodes(nodeToDistort2, 1) + s;
%         solidObject.meshObject.nodes(nodeToDistort2-1, 1) = solidObject.meshObject.nodes(nodeToDistort2-1, 1) + s/2;
%         solidObject.meshObject.nodes(nodeToDistort2+1, 1) = solidObject.meshObject.nodes(nodeToDistort2+1, 1) + s/2;
%         solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + lengthX/2;
%         solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2) + lengthY/2;

        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshRectangle(lengthX, lengthY, 2, 1, order, true);
        solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2) + lengthY/2;

        globalA = [0; 0];
        globalB = [1; 0];
        globalC = [2; 3];
        globalD = [-1; 2];

        R = [cos(phi), sin(phi); -sin(phi), cos(phi)];

        A = R*globalA;
        B = R*globalB;
        C = R*globalC;
        D = R*globalD;


%         [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity);

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
        
        % Nodal Forces
        nodalForceObject = nodalForceClass(dofObject);
        nodalForceObject.dimension = 2;                                % Dimension des zu betrachtenden Körpers angeben
        nodalForceObject.masterObject = solidObject;
        nodalForceObject.forceVector = [200; 0];                         % Komponenten des Kraftvektors anpassen
        nodalForceObject.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == lengthY);

        nodalForceObject2 = nodalForceClass(dofObject);
        nodalForceObject2.dimension = 2;                                % Dimension des zu betrachtenden Körpers angeben
        nodalForceObject2.masterObject = solidObject;
        nodalForceObject2.forceVector = -[200; 0];                         % Komponenten des Kraftvektors anpassen
        nodalForceObject2.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == 0);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        nodeToMeasure = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == 0);
        endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement;
        
    end
end