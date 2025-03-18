%% Start parameters
elementsToTest = {'BD', 'dispPetrovGalerkinBD'}; % {'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD'};
numberStepsMeshDistortion = 1;
maxDistortion = 5;

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
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'zero';
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-6;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        plateObject = plateClass(dofObject);
        switch elementsToTest{ii}
            case 'disp'
                plateObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'dispPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'displacementPetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'dispPetrovGalerkinBD'
                plateObject.materialObject.name = 'HookeBatheDvorkin';
                plateObject.elementDisplacementType = 'displacement';
                plateObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'BD'
                plateObject.materialObject.name = 'HookeBD';
            case 'eas'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'eas';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easBD'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'eas';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easPetrovGalerkin'
                plateObject.materialObject.name = 'Hooke';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'easPetrovGalerkinBD'
                plateObject.materialObject.name = 'HookeBD';
                plateObject.elementDisplacementType = 'easPetrovGalerkin';
                plateObject.mixedFEObject.condensation = true;
                plateObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        plateLengthX = 50;
        plateLengthY = 50;
        [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshPlateDistorted(plateLengthX, plateLengthY, s, order, serendipity);
%         [plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshPatchTestDistorted2D(plateLengthX, plateLengthY, order, serendipity);
        %plateObject.materialObject.name = 'Hooke';
        plateObject.materialObject.rho = 0;
        plateObject.materialObject.E = 1e4;
        plateObject.materialObject.nu = 0;
        plateObject.h = 1;
        plateObject.dimension = 2;
        plateObject.shapeFunctionObject.order = order;
        plateObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;

        % Dirichlet boundary
        boundary1 = dirichletClass(dofObject);
        boundary1.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        boundary1.nodalDof = 1;
        boundary1.masterObject = plateObject;
        
        boundary2 = dirichletClass(dofObject);
        boundary2.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
        boundary2.nodalDof = 2;
        boundary2.masterObject = plateObject;
        
%         boundary3 = dirichletClass(dofObject);
%         boundary3.nodeList = find(plateObject.meshObject.nodes(:,1) == 0);
%         boundary3.nodalDof = 3;
%         boundary3.masterObject = plateObject;
%         
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
%         neumannObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
%         neumannObject.masterObject = plateObject;
%         neumannObject.forceVector = [-1; 0; 0];                         % Komponenten des Kraftvektors anpassen
%         neumannObject.shapeFunctionObject.dimension = plateObject.dimension;
%         neumannObject.shapeFunctionObject.order = plateObject.shapeFunctionObject.order;
%         neumannObject.shapeFunctionObject.numberOfGausspoints = plateObject.shapeFunctionObject.numberOfGausspoints;  % Anzahl der Gaußpunkten auf dem Neumannrand anpassen
%         neumannObject.meshObject.edof = plateObject.meshObject.edof;
        
        % Nodal Forces
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = plateObject;
        nodalLoadObject.loadVector = -[0; 150; 0];                         % Komponenten des Kraftvektors anpassen
        nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX);
        
%         nodalLoadObject = nodalLoadClass(dofObject);
%         nodalLoadObject.masterObject = plateObject;
%         nodalLoadObject.loadVector = -[16.3666; 0; 0]/4;        % -[16.3527; 0; 0]/4                 % Komponenten des Kraftvektors anpassen
%         nodalLoadObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLengthX);
        
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        centralDisplacement = max(abs(plateObject.qN1(:, 1)));
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
ylim([min(min(resultArray)), max(max(resultArray))]);
legend(elementsToTest);

function [nodes, edof] = meshPlateDistorted(plateLengthX, plateLengthY, s, order, serendipity)
        [nodes,edof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
        nodeToDistort = find(nodes(:,1) == 0 & nodes(:,2) == 0);
%         nodes(nodeToDistort, :) = nodes(nodeToDistort, :) + s;
        nodes(nodeToDistort, :) = nodes(nodeToDistort, :) + s;
        if order == 2 && serendipity
            nodes(nodeToDistort-1, 1) = nodes(nodeToDistort-1, 1) + s/2;
            nodes(nodeToDistort+1, 1) = nodes(nodeToDistort+1, 1) + s/2;
            nodes(nodeToDistort+4, 1) = nodes(nodeToDistort+4, 1) + s/2;
            nodes(nodeToDistort-4, 1) = nodes(nodeToDistort-4, 1) + s/2;
        end
        nodes(:,1) = nodes(:,1) + plateLengthX/2;
        nodes(:,2) = nodes(:,2) + plateLengthY/2;
end