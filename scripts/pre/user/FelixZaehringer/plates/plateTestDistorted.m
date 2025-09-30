%% Start parameters
elementsToTest = {'easPetrovGalerkinFullIncompatible10Field', 'displacementBatheDvorkin'};%, 'displacementBatheDvorkin', 'displacementQ8', 'displacementPetrovGalerkinQ8', 'easAndelfingerRamm', 'easPetrovGalerkinAndelfingerRamm', 'easSimoRifai'
numberStepsMeshDistortion = 10;
maxDistortion = 10;
testBothDirections = false;

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
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress = struct('type', 'Cauchy', 'component', 13);
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-4;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        plateObject = plateClass(dofObject);
        numberOfGausspoints = @(order) (order + 1)^2;
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
            case 'selectiveReducedIntegrationQ8'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                plateObject.elementDisplacementType = 'selectiveReducedIntegration';
                serendipity = true;
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
            case 'displacementPetrovGalerkinAssumedStrain'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinAssumedStrain';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinBatheDvorkinMixed'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkinMixed';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinBatheDvorkinSkewSymbolic'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkinSkewSymbolic';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinBatheDvorkinSymbolic'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkinSymbolic';
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
                plateObject.elementNameAdditionalSpecification = 'BatheDvorkin1985';
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
            case 'easPetrovGalerkinIncompatible'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinIncompatible';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinFullIncompatible'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinFullIncompatible';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 6;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinFullIncompatible4DOF'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinFullIncompatible4DOF';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinFullIncompatible6DOF'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinFullIncompatible6DOF';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 6;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinFullIncompatible8DOF'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinFullIncompatible8DOF';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 8;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinFullIncompatible10Field'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinFullIncompatible10Field';
                plateObject.mixedFEObject.condensation = false;
                plateObject.mixedFEObject.typeShapeFunctionData = 10;
                order = 1;
                numberOfGausspoints = @(order) (order + 2)^2;
                serendipity = false;
            case 'easPetrovGalerkinAssumedStrain'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'AssumedStrain';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinAssumedStrain2iDOF'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'AssumedStrain2iDOF';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 2;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinAssumedStrain6iDOF'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'AssumedStrain6iDOF';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 6;
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
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'AndelfingerRamm';
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
            case 'easPetrovGalerkinBatheDvorkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin';
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

        plateLengthX = 50;
        plateLengthY = 50;
        [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshPlateDistorted(plateLengthX, plateLengthY, s, max(order), serendipity);

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof, bounEdof] = meshRectangle(plateLengthX, plateLengthY, 80, 80, order, serendipity);
%         plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
%         plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;

        %plateObject.materialObject.name = 'Hooke';
        plateObject.materialObject.rho = 0;
        plateObject.materialObject.E = 1e4;
        plateObject.materialObject.nu = 0.3;
        plateObject.h = 1;
        plateObject.dimension = 2;
        plateObject.shapeFunctionObject.order = order;
        plateObject.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints(order);

        % Dirichlet boundary
        boundary1 = dirichletClass(dofObject);
        boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0);
        boundary1.nodalDof = 1;
        boundary1.masterObject = plateObject;
        
        boundary2 = dirichletClass(dofObject);
        boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,1) == plateLengthX);
        boundary2.nodalDof = 2;
        boundary2.masterObject = plateObject;
        
        boundary3 = dirichletClass(dofObject);
        boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0 | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
        boundary3.nodalDof = 3;
        boundary3.masterObject = plateObject;
        
%         boundary1 = dirichletClass(dofObject);
%         boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
%         boundary1.nodalDof = 1;
%         boundary1.masterObject = plateObject;
%         
%         boundary2 = dirichletClass(dofObject);
%         boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
%         boundary2.nodalDof = 2;
%         boundary2.masterObject = plateObject;
%         
%         boundary3 = dirichletClass(dofObject);
%         boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
%         boundary3.nodalDof = 3;
%         boundary3.masterObject = plateObject;
        
        % % Neumann boundary
%         neumannObject = neumannClass(dofObject);
%         neumannObject.masterObject = plateObject;
%         neumannObject.loadGeometry = 'area';
%         neumannObject.loadVector = [-4; 10; 0];                         % Komponenten des Kraftvektors anpassen
%         neumannObject.meshObject.edof = plateObject.meshObject.edof;
%         
        % Nodal Forces
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = plateObject;
        nodalLoadObject.loadVector = -[16.3527; 0; 0]/4;        % -[16.3527; 0; 0]/4       -[16.3666; 0; 0]/4          % Komponenten des Kraftvektors anpassen
        nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX & plateObject.meshObject.nodes(:,2) == plateLengthY);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        centralDisplacement = max(abs(plateObject.qN1(nodalLoadObject.nodeList, 1)));
        disp(['Central displacement: ', num2str(centralDisplacement)]);
        resultArray(jj, ii) = centralDisplacement;
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
% ylim([min(min(resultArray)), max(max(resultArray))]);
% ylim([0.0, 0.9]);
xlabel('s');
ylabel('w_{max}');
legend(elementsToTest);

function [nodes, edof] = meshPlateDistorted(plateLengthX, plateLengthY, s, order, serendipity)
        [nodes, edof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
        nodeToDistort = find(nodes(:,1) == 0 & nodes(:,2) == 0);
        nodes(nodeToDistort, :) = nodes(nodeToDistort, :) + s;
        nodes = adaptMeshToCornerNodes(nodes, edof, 2);
%         nodes = adaptMeshToCornerNodes(nodes, edof, 2);
        nodes(:,1) = nodes(:,1) + plateLengthX/2;
        nodes(:,2) = nodes(:,2) + plateLengthY/2;
end

function [nodes, edof] = meshPlateDistortedAlternative(plateLengthX, plateLengthY, s, order, serendipity)
        [nodes, edof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
        nodesToDistort = find(nodes(:,1) == 0);
        nodes(nodesToDistort, 1) = nodes(nodesToDistort, 1) + s;
%         nodes = adaptMeshToCornerNodes(nodes, edof, 2);
        nodes(:,1) = nodes(:,1) + plateLengthX/2;
        nodes(:,2) = nodes(:,2) + plateLengthY/2;
end