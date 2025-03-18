%% Start parameters
elementsToTest = {'eas'}; % {'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD'};
lengthX = 50;
lengthY = 5;
numberStepsMeshDistortion = 10;
maxDistortion = min(lengthX/2, lengthY/2)*0.8;
testBothDirections = true;

%% Script
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
                elementDisplacementType = 'easPetrovGalerkin';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end
        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshPlateDistorted(lengthX, lengthY, s, order, serendipity);
        solidObject.elementDisplacementType = elementDisplacementType;
        solidObject.materialObject.name = materialName;

        %plateObject.materialObject.name = 'Hooke';
        solidObject.materialObject.rho = 0;
        solidObject.materialObject.E = 1e4;
        solidObject.materialObject.nu = 0.3;
        solidObject.dimension = 2;
        solidObject.shapeFunctionObject.order = order;
        solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
        dirichletBoundary.nodalDof = 1:2;
        dirichletBoundary.masterObject = solidObject;
        
        % Nodal Forces
        nodalLoadObject = nodalLoadClass(dofObject);
        nodalLoadObject.masterObject = solidObject;
        nodalLoadObject.loadVector = [0; 10];                         % Komponenten des Kraftvektors anpassen
        nodalLoadObject.nodeList = find(solidObject.meshObject.nodes(:,1) == lengthX);
        
        
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

function [nodes, edof] = meshPlateDistorted(plateLengthX, plateLengthY, s, order, serendipity)
        [nodes, edof] = meshRectangle(plateLengthX, plateLengthY, 2, 2, order, serendipity);
        nodeToDistort = find(nodes(:,1) == 0 & nodes(:,2) == 0);
        nodes(nodeToDistort, :) = nodes(nodeToDistort, :) + s;
        nodes = adaptMeshToCornerNodes(nodes, edof, 2);
        nodes(:,1) = nodes(:,1) + plateLengthX/2;
        nodes(:,2) = nodes(:,2) + plateLengthY/2;
end