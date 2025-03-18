%% Start parameters
elementsToTest = {'displacementBatheDvorkin', 'displacementPetrovGalerkinBatheDvorkin', 'displacementQ8', 'displacementPetrovGalerkinQ8'};
numberStepsMeshDistortion = 20;
maxDistortion = 4;
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
        setupObject.saveObject.fileName = 'plateSimple';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.usePreconditioning = false;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        plateLengthX = 10;
        plateLengthY = 10;
        plateObject = plateClass(dofObject);

        switch elementsToTest{ii}
            case 'displacementQ8'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'displacementQ4'
                plateObject.materialObject.name = 'Hooke';
                order = 1;
                serendipity = false;
            case 'selectiveReducedIntegration'
                plateObject.materialObject.name = 'Hooke';
                order = 1;
                plateObject.elementDisplacementType = 'selectiveReducedIntegration';
                serendipity = false;
                plateObject.mixedFEObject.condensation = true;
            case 'selectiveReducedIntegrationPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                order = 1;
                plateObject.elementDisplacementType = 'selectiveReducedIntegration';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                serendipity = false;
                plateObject.mixedFEObject.shapeFunctionObject.numberOfGausspoints = 4;
                plateObject.mixedFEObject.condensation = true;
            case 'displacementPetrovGalerkinQ4'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'displacementQ9'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = false;
            case 'displacementQ16'
                plateObject.materialObject.name = 'Hooke';
                order = 3;
                serendipity = false;
            case 'displacementPetrovGalerkinQ8'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'displacementPetrovGalerkinQ9'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = false;
            case 'displacementPetrovGalerkinBatheDvorkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin';
                order = 1;
                serendipity = false;
            case 'displacementPartialPetrovGalerkinBatheDvorkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PartialPetrovGalerkinBatheDvorkin';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinWOnlyBatheDvorkin'
                plateObject.materialObject.name = 'WOnlyHookeBD';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'displacementBatheDvorkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
                order = 1;
                serendipity = false;
            case 'displacementHughesTezduyar'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementNameAdditionalSpecification = 'HughesTezduyar';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinHughesTezduyar'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinHughesTezduyar';
                order = 1;
                serendipity = false;
            case 'easSimoRifai'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'mixed'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'mixed';
                plateObject.elementNameAdditionalSpecification = '';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'mixedQ8'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'mixed';
                plateObject.elementNameAdditionalSpecification = '';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 10;
                order = 2;
                serendipity = true;
            case 'mixedPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'mixed';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easSimoRifaiE8'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 8;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinSimoRifai'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinSimoRifaiE6'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 6;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinSimoRifaiPurelyLagrangeE6'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'SimoRifaiPurelyLagrange';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 6;
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinSimoRifaiBbarE6'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacementPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'SimoRifaiBbar';
                order = 1;
                serendipity = false;
            case 'easAndelfingerRamm'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'eas';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinAndelfingerRamm'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinAndelfingerRammSelectiveReducedIntegration'
                plateObject.materialObject.name = 'Hooke';
                order = 1;
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinSelectiveReducedIntegration';
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                serendipity = false;
                plateObject.mixedFEObject.condensation = true;
            case 'mixedDisplacementHooke'
                plateObject.materialObject.name = 'Hooke';
                numberOfNodes = [9, 4, 4];
                order = 2;
                plateObject.shapeFunctionObject.numberOfNodes = numberOfNodes;
                plateObject.elementNameAdditionalSpecification = 'Mixed';
                serendipity = false;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, 2, 1, order, serendipity);
        nodeToDistort = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == -plateLengthY/2);
        nodeToDistort2 = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == +plateLengthY);
        plateObject.meshObject.nodes(nodeToDistort, 1) = plateObject.meshObject.nodes(nodeToDistort, 1) + s;
        plateObject.meshObject.nodes(nodeToDistort2, 1) = plateObject.meshObject.nodes(nodeToDistort2, 1) - s;
        if order==2 && serendipity
            plateObject.meshObject.nodes(nodeToDistort-1, 1) = plateObject.meshObject.nodes(nodeToDistort-1, 1) + s/2;
            plateObject.meshObject.nodes(nodeToDistort+1, 1) = plateObject.meshObject.nodes(nodeToDistort+1, 1) + s/2;
            plateObject.meshObject.nodes(nodeToDistort2-1, 1) = plateObject.meshObject.nodes(nodeToDistort2-1, 1) - s/2;
            plateObject.meshObject.nodes(nodeToDistort2+1, 1) = plateObject.meshObject.nodes(nodeToDistort2+1, 1) - s/2;
        end
        plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
        plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshPatchTestDistorted2D(plateLengthX, plateLengthY, order, serendipity);

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
%         nodeToDistort = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == 0);
%         plateObject.meshObject.nodes(nodeToDistort, 1) = plateObject.meshObject.nodes(nodeToDistort, 1) + s;
%         if order==2 && serendipity
%             plateObject.meshObject.nodes(nodeToDistort-1, 1) = plateObject.meshObject.nodes(nodeToDistort-1, 1) + s/2;
%             plateObject.meshObject.nodes(nodeToDistort+1, 1) = plateObject.meshObject.nodes(nodeToDistort+1, 1) + s/2;
%             plateObject.meshObject.nodes(nodeToDistort-4, 1) = plateObject.meshObject.nodes(nodeToDistort-4, 1) + s/2;
%             plateObject.meshObject.nodes(nodeToDistort+4, 1) = plateObject.meshObject.nodes(nodeToDistort+4, 1) + s/2;
%         end
%         plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
%         plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, 1, 2, order, serendipity);
%         nodeToDistort1 = find(plateObject.meshObject.nodes(:,1) == -plateLengthX/2 & plateObject.meshObject.nodes(:,2) == 0);
%         nodeToDistort2 = find(plateObject.meshObject.nodes(:,1) == plateLengthX/2 & plateObject.meshObject.nodes(:,2) == 0);
%         plateObject.meshObject.nodes(nodeToDistort1, 2) = plateObject.meshObject.nodes(nodeToDistort1, 2) + s;
%         plateObject.meshObject.nodes(nodeToDistort2, 2) = plateObject.meshObject.nodes(nodeToDistort2, 2) - s;
%         if order==2 && serendipity
%             plateObject.meshObject.nodes(nodeToDistort1-2, 2) = plateObject.meshObject.nodes(nodeToDistort1-2, 2) + s/2;
%             plateObject.meshObject.nodes(nodeToDistort1+3, 2) = plateObject.meshObject.nodes(nodeToDistort1+3, 2) + s/2;
%             plateObject.meshObject.nodes(nodeToDistort2-3, 2) = plateObject.meshObject.nodes(nodeToDistort2-3, 2) - s/2;
%             plateObject.meshObject.nodes(nodeToDistort2+2, 2) = plateObject.meshObject.nodes(nodeToDistort2+2, 2) - s/2;
%         end
%         plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
%         plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;
        

        %plateObject.materialObject.name = 'Hooke';
        plateObject.materialObject.rho = 0;
        plateObject.materialObject.E = 1e5;
        plateObject.materialObject.nu = 0;
        plateObject.h = 1;
        plateObject.dimension = 2;
        plateObject.shapeFunctionObject.order = order;
        plateObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        boundary1 = dirichletClass(dofObject);
        boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        boundary1.nodalDof = 1:3;
        boundary1.masterObject = plateObject;
        
        % boundary1 = dirichletClass(dofObject);
        % boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        % boundary1.nodalDof = 1;
        % boundary1.masterObject = plateObject;
        % 
        % boundary2 = dirichletClass(dofObject);
        % boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        % boundary2.nodalDof = 2;
        % boundary2.masterObject = plateObject;
        % 
        % boundary3 = dirichletClass(dofObject);
        % boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        % boundary3.nodalDof = 3;
        % boundary3.masterObject = plateObject;
        
        % % Neumann boundary
%         neumannObject = neumannClass(dofObject);
%         neumannObject.masterObject = plateObject;
%         neumannObject.loadGeometry = 'area';
%         neumannObject.loadVector = [-1; 0; 0];                         % Komponenten des Kraftvektors anpassen
%         neumannObject.meshObject.edof = plateObject.meshObject.edof;
        
        % Nodal Forces
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = plateObject;
        nodalLoadObject.loadVector = [1; 0; 0];
        nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX & plateObject.meshObject.nodes(:,2) == plateLengthY);
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        endDisplacement = max(abs(plateObject.qN1(:, 1)));
        disp(['Central displacement: ', num2str(endDisplacement)]);
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