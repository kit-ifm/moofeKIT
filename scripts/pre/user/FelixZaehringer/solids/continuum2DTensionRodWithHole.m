% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'dispQ4'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
numberOfRefinements = 0;

%% Script
resultArray = zeros(numberOfRefinements+1, size(elementsToTest, 2));
hArray = zeros(numberOfRefinements+1, 1);
numberOfDofArray = zeros(numberOfRefinements+1, 1);
numberOfElementsArray = zeros(numberOfRefinements+1, 1);
for ii = 1:size(elementsToTest, 2)
    for jj = 1:numberOfRefinements+1
        radius = 2;
        lengthX = 6;
        lengthY = 12;

        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp(['Current refinement: ', num2str(jj)]);
        disp('=======================================');

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = true;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 22;
        setupObject.plotObject.view = 2;
        setupObject.usePreconditioning = false;
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';
        
        dofObject = dofClass;   % required object for dof and object handling

        %% continuum Objects
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
                solidObject.elementDisplacementType = 'easPetrovGalerkin';
                solidObject.mixedFEObject.condensation = false;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end

        [solidObject.meshObject.nodes,solidObject.meshObject.edof, bounEdof] = meshTensionRodWithHole(radius, lengthX, lengthY, jj*10, order, serendipity);

        numberOfDofArray(jj) = length(unique(solidObject.meshObject.edof));
        numberOfElementsArray (jj) = size(solidObject.meshObject.edof, 1);

        %plateObject.materialObject.name = 'Hooke';
        solidObject.materialObject.rho = 0;
        solidObject.materialObject.E = 210000;
        solidObject.materialObject.nu = 0.3;
        solidObject.dimension = 2;
        solidObject.shapeFunctionObject.order = order;
        solidObject.shapeFunctionObject.numberOfGausspoints = (order + 1)^2;
        
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
        dirichletBoundary.nodalDof = 1;
        dirichletBoundary.masterObject = solidObject;

        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
        dirichletBoundary2.nodalDof = 2;
        dirichletBoundary2.masterObject = solidObject;
        
        % Neumann boundary
        fy = 1;
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) [0; fy]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof;
        
        %% solver
        dofObject = runNewton(setupObject, dofObject);
        evalPointX = radius;
        evalPointY = 0;
        nodeToMeasure = find(abs(solidObject.meshObject.nodes(:,1) - evalPointX) < 1e-8 & abs(solidObject.meshObject.nodes(:,2) - evalPointY) < 1e-8);
        stress = dofObject.listContinuumObjects{1}.plotData.nodalColorData(nodeToMeasure);
        disp(['Current Stress: ', num2str(stress)]);
        nennspannung = fy*lengthX/(lengthX - radius);
        disp(['Ratio: ', num2str(stress/nennspannung)]);
    end
end

% evalStart = 1;
% xVals = numberOfDofArray;
% figure;
% for ii=1:size(elementsToTest, 2)
%     fit=polyfit(log(xVals(evalStart:end)),log(resultArray(evalStart:end, ii)),1);
%     x1=min(xVals(evalStart:end)):0.1:max(xVals(evalStart:end));
%     y1=exp(fit(1)*log(x1)+fit(2)); 
%     loglog(xVals(evalStart:end), resultArray(evalStart:end, ii), 'o', x1, y1, '-');
%     loglog(x1, y1, '-');
%     hold on;
%     disp(['Steigung: ', num2str(fit(1))])
% end
% xlabel('s');
% ylabel('w_{max}');
% legend(elementsToTest);
