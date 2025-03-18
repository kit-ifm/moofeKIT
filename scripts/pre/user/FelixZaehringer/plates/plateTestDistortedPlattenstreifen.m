%% Start parameters
elementsToTest = {'BD'}; % {'BD', 'dispPetrovGalerkinBD', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkin'};
numberStepsMeshDistortion = 1;
maxDistortion = 2;
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
        setupObject.saveObject.saveData = true;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 23;
        setupObject.usePreconditioning = false;
%         setupObject.newton.tolerance = 1e-7;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        plateLengthX = 10;
        plateLengthY = 10;
        plateObject = plateClass(dofObject);

        initialX = 0;

        switch elementsToTest{ii}
            case 'dispQ8'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'dispQ4'
                plateObject.materialObject.name = 'Hooke';
                order = 1;
                serendipity = false;
            case 'dispQ9'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = false;
            case 'dispPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'dispPetrovGalerkinQ4'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'dispPetrovGalerkinBD'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin';
                order = 1;
                serendipity = false;
            case 'BD'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementNameAdditionalSpecification = 'BatheDvorkin';
                order = 1;
                serendipity = false;
            case 'eas'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easBD'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'eas';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.elementNameAdditionalSpecification = 'SimoRifai';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            case 'easPetrovGalerkinBD'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
                order = 1;
                serendipity = false;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof, bounEdof] = meshRectangle(plateLengthX, plateLengthY, 2, 10, order, serendipity);

        [plateObject.meshObject.nodes, plateObject.meshObject.edof, bounEdof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
        nodeToDistort = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == 0);
%         nodeToDistort2 = find(plateObject.meshObject.nodes(:,1) == 0 & plateObject.meshObject.nodes(:,2) == +plateLengthY/2);
        plateObject.meshObject.nodes(nodeToDistort, 1) = plateObject.meshObject.nodes(nodeToDistort, 1) + s;
        plateObject.meshObject.nodes(nodeToDistort, 2) = plateObject.meshObject.nodes(nodeToDistort, 2) + s;
%         plateObject.meshObject.nodes = adaptMeshToCornerNodes(plateObject.meshObject.nodes, plateObject.meshObject.edof, 2);
        plateObject.meshObject.nodes(:,1) = plateObject.meshObject.nodes(:,1) + plateLengthX/2;
        plateObject.meshObject.nodes(:,2) = plateObject.meshObject.nodes(:,2) + plateLengthY/2;

%         [plateObject.meshObject.nodes,plateObject.meshObject.edof, bounEdof] = meshPatchTestDistorted2D(plateLengthX, plateLengthY, order, serendipity);


%         [plateObject.meshObject.nodes,plateObject.meshObject.edof, bounEdof] = meshPlateStripDistorted(plateLengthX, plateLengthY, 2, 20, order, serendipity, s);

        plateObject.meshObject.nodes(:, 1) = plateObject.meshObject.nodes(:, 1) + initialX;

        %plateObject.materialObject.name = 'Hooke';
        plateObject.materialObject.rho = 0;
        plateObject.materialObject.E = 1200;
        plateObject.materialObject.nu = 0;
        if plateObject.materialObject.nu ~= 0
            plateObject.materialObject.E = 1200*(1-plateObject.materialObject.nu^2);
        end
        plateObject.h = 1;
        plateObject.dimension = 2;
        plateObject.shapeFunctionObject.order = order;
        plateObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        boundary1 = dirichletClass(dofObject);
        boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == initialX);
        boundary1.nodalDof = 1:3;
        boundary1.masterObject = plateObject;
        
        % boundary1 = dirichletClass(dofObject);
        % boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        % boundary1.nodalDof = 1;
        % boundary1.masterObject = plateObject;
        % 
%         boundary2 = dirichletClass(dofObject);
%         boundary2.nodeList = find(plateObject.meshObject.nodes(:,2) == plateLengthY/2 | plateObject.meshObject.nodes(:,2) == -plateLengthY/2);
%         boundary2.nodalDof = 2;
%         boundary2.masterObject = plateObject;
        % 
        % boundary3 = dirichletClass(dofObject);
        % boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        % boundary3.nodalDof = 3;
        % boundary3.masterObject = plateObject;
        
        % Neumann boundary
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = plateObject;
        neumannObject.loadVector = [2; 0; 0];                         % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX2;

%         neumannObject = neumannClass(dofObject);
%         neumannObject.loadGeometry = 'line';
%         neumannObject.masterObject = plateObject;
%         neumannObject.loadVector = [0; 2; 0];                         % Komponenten des Kraftvektors anpassen
%         neumannObject.meshObject.edof = bounEdof.SX2;

%         neumannObject = neumannClass(dofObject);
%         neumannObject.loadGeometry = 'area';
%         neumannObject.masterObject = plateObject;
%         neumannObject.loadVector = [2; 0; 0];                         % Komponenten des Kraftvektors anpassen
%         neumannObject.meshObject.edof = plateObject.meshObject.edof;
        
%         % Nodal Forces
%         nodalForceObject = nodalForceClass(dofObject);
%         nodalForceObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
%         nodalForceObject.masterObject = plateObject;
%         nodalForceObject.forceVector = [0; 20; 0];                         % Komponenten des Kraftvektors anpassen
% %         nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY));
%         nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        centralNode = find(plateObject.meshObject.nodes(:,1) == plateLengthX+initialX & plateObject.meshObject.nodes(:,2) == plateLengthY/2);
        endDisplacement = plateObject.qN1(centralNode, 1);
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
% ylim([min(min(resultArray)), max(max(resultArray))]);
% ylim([0.0, 0.9]);
xlabel('s');
ylabel('w_{max}');
legend(elementsToTest);

function [nodes, edof, bounEdof] = meshPlateStripDistorted(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity, s)
    assert(mod(numberOfElementsX, 2) == 0, 'Number of elements in X direction must be even!');
    assert(mod(order, 2) == 0 | (order == 1 & mod(numberOfElementsY, 2) == 0), 'Number of elements in Y direction must be even for uneven order!');
    elementLengthY = lengthY/numberOfElementsY;
    cornerNodesY = -lengthY/2:elementLengthY:lengthY/2;
    [nodes, edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
    nodesToDistort = find(nodes(:, 1) == 0 & ismember(nodes(:, 2), cornerNodesY));
    for ii=1:length(nodesToDistort)
        currentNode = nodesToDistort(ii);
        if mod(ii, 2) ~= 0
            nodes(currentNode, 1) = nodes(currentNode, 1) + s;
        else
            nodes(currentNode, 1) = nodes(currentNode, 1) - s;
%             nodes(currentNode, 2) = nodes(currentNode, 2) - s;
        end
    end
    nodes = adaptMeshToCornerNodes(nodes, edof, 2);
    nodes(:,1) = nodes(:,1) + lengthX/2;
    nodes(:,2) = nodes(:,2) + lengthY/2;
end