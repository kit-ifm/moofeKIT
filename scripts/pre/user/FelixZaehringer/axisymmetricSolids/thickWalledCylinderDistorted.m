%% Start parameters
elementsToTest = {'displacementQ8', 'displacementPetrovGalerkinQ8'}; % 'displacementPetrovGalerkinQ8'
numberStepsMeshDistortion = 100;
maxDistortion = 0.9;
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
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 22;
        setupObject.plotObject.view = 2;
        setupObject.integrator = 'Endpoint';
        setupObject.newton.tolerance = 1e-4;
        
        dofObject = dofClass;   % required object for dof and object handling
        
        %% continuum Objects
        axisymmetricSolidObject = axisymmetricSolidClass(dofObject);
        switch elementsToTest{ii}
            case 'displacementQ8'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = true;
            case 'displacementQ4'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 1;
                serendipity = false;
            case 'displacementPetrovGalerkinQ4'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 1;
                serendipity = false;
            case 'displacementQ9'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                order = 2;
                serendipity = false;
            case 'displacementPetrovGalerkinQ8'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'displacement';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = true;
            case 'displacementPetrovGalerkinQ9'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'displacement';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'PetrovGalerkin';
                order = 2;
                serendipity = false;
            case 'easSimoRifai'
                axisymmetricSolidObject.materialObject.name = 'Hooke';
                axisymmetricSolidObject.elementDisplacementType = 'eas';
                axisymmetricSolidObject.elementNameAdditionalSpecification = 'SimoRifai';
                axisymmetricSolidObject.mixedFEObject.condensation = true;
                axisymmetricSolidObject.mixedFEObject.typeShapeFunctionData = 5;
                order = 1;
                serendipity = false;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        lengthX = 4;
        lengthY = 4;
        numberOfElementsX = 2;
        numberOfElementsY = 2;
        innerRadius = 1;
        outerRadius = innerRadius + lengthX;
        [axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity);
        nodeToDistort1 = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
        nodeToDistort2 = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
        nodeToDistort = find(axisymmetricSolidObject.meshObject.nodes(:, 1) == 0 & axisymmetricSolidObject.meshObject.nodes(:, 2) == 0);
        axisymmetricSolidObject.meshObject.nodes(:, 1) = axisymmetricSolidObject.meshObject.nodes(:, 1) + lengthX/2 + innerRadius;

%         axisymmetricSolidObject.meshObject.nodes(nodeToDistort1, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort1, 1) + s;
%         axisymmetricSolidObject.meshObject.nodes(nodeToDistort2, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort2, 1) - s;
        axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) + s;
%         axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 1) + s;
        axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) = axisymmetricSolidObject.meshObject.nodes(nodeToDistort, 2) + s;
        axisymmetricSolidObject.meshObject.nodes = adaptMeshToCornerNodes(axisymmetricSolidObject.meshObject.nodes, axisymmetricSolidObject.meshObject.edof, 2);

        axisymmetricSolidObject.materialObject.rho = 0;
        axisymmetricSolidObject.materialObject.E = 1e7;
        axisymmetricSolidObject.materialObject.nu = 0.3;
        axisymmetricSolidObject.dimension = 2;
        axisymmetricSolidObject.shapeFunctionObject.order = order;
        axisymmetricSolidObject.shapeFunctionObject.numberOfGausspoints = (max(order) + 1)^2;

        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == -lengthY/2);
        dirichletBoundary.nodalDof = 2;
        dirichletBoundary.timeFunction = @(t, XYZ) XYZ;
        dirichletBoundary.masterObject = axisymmetricSolidObject;
        
        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(axisymmetricSolidObject.meshObject.nodes(:, 2) == lengthY/2);
        dirichletBoundary2.nodalDof = 2;
        dirichletBoundary2.timeFunction = @(t, XYZ) lengthY/2;
        dirichletBoundary2.masterObject = axisymmetricSolidObject;


        % Neumann boundary
        evalPointX = innerRadius;
        evalPointY = -lengthY/2;
        p = axisymmetricSolidObject.materialObject.E*(outerRadius^2-innerRadius^2)/(1+axisymmetricSolidObject.materialObject.nu)/innerRadius^2/(outerRadius^2/evalPointX+(1-2*axisymmetricSolidObject.materialObject.nu)*evalPointX);
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = axisymmetricSolidObject;
        neumannObject.loadVector = [p; 0]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX1;
        
        
        %% solver
        dofObject = runNewton(setupObject, dofObject);
        nodeToMeasure = find(abs(axisymmetricSolidObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(axisymmetricSolidObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8);
        endDisplacement = axisymmetricSolidObject.qN1(nodeToMeasure, 1) - axisymmetricSolidObject.qN(nodeToMeasure, 1);
        disp(['End displacement: ', num2str(endDisplacement)]);
        exactDisplacement = (1+axisymmetricSolidObject.materialObject.nu)*p*innerRadius^2/(axisymmetricSolidObject.materialObject.E*(outerRadius^2-innerRadius^2))*(outerRadius^2/evalPointX+(1-2*axisymmetricSolidObject.materialObject.nu)*evalPointX);
        disp(['Ratio FEM/Exact: ', num2str(endDisplacement/exactDisplacement)]);
        resultArray(jj, ii) = endDisplacement/exactDisplacement;
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
ylabel('u(r_i)');
legend(elementsToTest);
