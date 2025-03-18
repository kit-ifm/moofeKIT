%% Start parameters
elementsToTest = {'easPetrovGalerkin', 'dispQ8'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD'};
lengthX = 10;
lengthY = 2;
numberStepsMeshDistortion = 50;
maxDistortion = 4.8;
testBothDirections = true;
%% Script
totalNumberOfSteps = numberStepsMeshDistortion+1;
if testBothDirections
    totalNumberOfSteps = 2*numberStepsMeshDistortion+1;
end
resultArray = zeros(totalNumberOfSteps, size(elementsToTest, 2));
for ii = 1:size(elementsToTest, 2)
    for jj = 1:totalNumberOfSteps
        s = (jj - 1)*maxDistortion/numberStepsMeshDistortion;
        if testBothDirections
            s = (jj - 1 - numberStepsMeshDistortion)*maxDistortion/numberStepsMeshDistortion;
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
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.usePreconditioning = false;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        solidObject = solidClass(dofObject);
        order = 1;
        serendipity = false;
        switch elementsToTest{ii}
            case 'disp'
                materialName = 'HookeESZ';
                elementDisplacementType = 'displacement';
            case 'dispQ9'
                materialName = 'HookeESZ';
                elementDisplacementType = 'displacement';
                order = 2;
            case 'dispQ8'
                materialName = 'HookeESZ';
                elementDisplacementType = 'displacement';
                order = 2;
                serendipity = true;
            case 'dispPetrovGalerkin'
                materialName = 'HookeESZ';
                elementDisplacementType = 'displacementPetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'eas'
                materialName = 'HookeESZ';
                elementDisplacementType = 'eas';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easPetrovGalerkin'
                materialName = 'HookeESZ';
                elementDisplacementType = 'eas';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end
        solidObject.elementDisplacementType = elementDisplacementType;
        solidObject.materialObject.name = materialName;

        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshRectangle(lengthX, lengthY, 2, 1, order, serendipity);
        nodeToDistort1 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY/2);
        nodeToDistort2 = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == -lengthY/2);
        solidObject.meshObject.nodes(nodeToDistort1, 1) = solidObject.meshObject.nodes(nodeToDistort1, 1) + s;
        solidObject.meshObject.nodes(nodeToDistort2, 1) = solidObject.meshObject.nodes(nodeToDistort2, 1) - s;
        if order == 2 && serendipity
        solidObject.meshObject.nodes(nodeToDistort1-1, 1) = solidObject.meshObject.nodes(nodeToDistort1-1, 1) + s/2;
        solidObject.meshObject.nodes(nodeToDistort1+1, 1) = solidObject.meshObject.nodes(nodeToDistort1+1, 1) + s/2;
        solidObject.meshObject.nodes(nodeToDistort2-1, 1) = solidObject.meshObject.nodes(nodeToDistort2-1, 1) - s/2;
        solidObject.meshObject.nodes(nodeToDistort2+1, 1) = solidObject.meshObject.nodes(nodeToDistort2+1, 1) - s/2;
        end
        solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2) + lengthY/2;

%         [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity);

        solidObject.materialObject.rho = 0;
        solidObject.materialObject.E = 1500;
        solidObject.materialObject.nu = 0.25;
        solidObject.dimension = 2;
        solidObject.shapeFunctionObject.order = order;
        solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary.nodalDof = 1:2;
        dirichletBoundary.masterObject = solidObject;

        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,1) == 0 & solidObject.meshObject.nodes(:,2) == lengthY);
        dirichletBoundary.nodalDof = 1;
        dirichletBoundary.masterObject = solidObject;
        
        % Nodal Forces
%         nodalForceObject = nodalForceClass(dofObject);
%         nodalForceObject.dimension = 2;                                % Dimension des zu betrachtenden Körpers angeben
%         nodalForceObject.masterObject = solidObject;
%         nodalForceObject.forceVector = -[10; 0];                         % Komponenten des Kraftvektors anpassen
%         nodalForceObject.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == lengthY/2);
% 
%         nodalForceObject2 = nodalForceClass(dofObject);
%         nodalForceObject2.dimension = 2;                                % Dimension des zu betrachtenden Körpers angeben
%         nodalForceObject2.masterObject = solidObject;
%         nodalForceObject2.forceVector = [10; 0];                         % Komponenten des Kraftvektors anpassen
%         nodalForceObject2.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == -lengthY/2);

        % Test Vertical Load
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = solidObject;
        nodalLoadObject.loadVector = [-10; 0];                         % Komponenten des Kraftvektors anpassen
        nodalLoadObject.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == lengthY);

        nodalLoadObject2 = nodalLoadClass(dofObject);
        nodalLoadObject2.masterObject = solidObject;
        nodalLoadObject2.loadVector = [10; 0];                         % Komponenten des Kraftvektors anpassen
        nodalLoadObject2.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX & solidObject.meshObject.nodes(:,2) == 0);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        endDisplacement = solidObject.qN1(end, 2) - solidObject.qN(end, 2);
        disp(['End displacement: ', num2str(endDisplacement)]);
        resultArray(jj, ii) = endDisplacement;
        
    end
end
figure;
if testBothDirections
    xVals = -maxDistortion:maxDistortion/numberStepsMeshDistortion:maxDistortion;
else
    xVals = 0:maxDistortion/numberStepsMeshDistortion:maxDistortion;
end

for ii=1:size(elementsToTest, 2)
        plot(xVals, resultArray(:, ii)');
    hold on;
end
xlim([min(xVals), max(xVals)]);
ylim([min(min(resultArray)), max(max(resultArray))]);
legend(elementsToTest);