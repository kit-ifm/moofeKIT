% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'huangEtAl2020'}; % {'disp', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
nuValuesToTest = [0.3, 0.49, 0.499, 0.4999, 0.49999];
%% Script
for ii = 1:size(elementsToTest, 2)
    for jj=1:length(nuValuesToTest)
        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp(['Current nu: ', num2str(nuValuesToTest(jj))]);
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
        setupObject.newton.tolerance = 1e-4;
        setupObject.integrator = 'Endpoint';
    
        dofObject = dofClass; % required object for dof and object handling
    
        %% continuum Objects
        solidObject = solidClass(dofObject);
    
        numberOfGausspoints = @(order) (order + 1)^2;

        switch elementsToTest{ii}
            case 'disp'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeEVZ';
            case 'dispQ4'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
            case 'dispQ9'
                order = 2;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
            case 'eas'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'eas';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkin'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkin';
            case 'dispPetrovGalerkinXie'
                order = 2;
                serendipity = true;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'displacementPetrovGalerkinXie';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'dispPetrovGalerkinCen'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'displacement';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinCenEtAl2015';
            case 'easPetrovGalerkin'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinPfefferkornBetsch2021';
                solidObject.elementDisplacementType = 'eas';
                solidObject.mixedFEObject.condensation = false;
                solidObject.mixedFEObject.typeShapeFunctionData = 4;
            case 'incompressibleHughes'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'incompressibleHughes';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 0;
            case 'incompressibleHerrmann'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'incompressibleHerrmann';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 0;
            case 'incompressibleSimoTaylorPistor'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeEVZ';
                solidObject.elementDisplacementType = 'incompressibleSimoTaylorPistor';
                solidObject.mixedFEObject.condensation = true;
                solidObject.mixedFEObject.typeShapeFunctionData = 0;
            case 'huangEtAl2020'
                order = 1;
                serendipity = false;
                solidObject.materialObject.name = 'HookeESZ';
                solidObject.elementDisplacementType = 'incompatibleModes';
                solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinHuangEtAl2020';
                solidObject.mixedFEObject.typeShapeFunctionData = 10;
                numberOfGausspoints = @(order) (order + 2)^2;
            otherwise
                warning(['Element ', elementsToTest{ii}, ' not implemented!'])
        end
    
        [solidObject.meshObject.nodes, solidObject.meshObject.edof, bounEdof] = meshThickWalledCylinder(order, serendipity);
    
        %plateObject.materialObject.name = 'Hooke';
        solidObject.materialObject.rho = 0;
        solidObject.materialObject.E = 1e3;
        solidObject.materialObject.nu = nuValuesToTest(jj);
        solidObject.dimension = 2;
        solidObject.shapeFunctionObject.order = order;
        solidObject.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints(order);
    
        % Dirichlet boundary
        dirichletBoundary = dirichletClass(dofObject);
        dirichletBoundary.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
        dirichletBoundary.nodalDof = 1;
        dirichletBoundary.masterObject = solidObject;
    
        dirichletBoundary2 = dirichletClass(dofObject);
        dirichletBoundary2.nodeList = find(solidObject.meshObject.nodes(:, 2) == 0);
        dirichletBoundary2.nodalDof = 2;
        dirichletBoundary2.masterObject = solidObject;
    
        % Neumann boundary
        p = 1;
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) p*[XYZ(1)/3; XYZ(2)/3]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof;
    
        %% solver
        dofObject = runNewton(setupObject, dofObject);
        evalPointX = 3;
        evalPointY = 0;
        nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-evalPointX) < 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-evalPointY) < 1e-8);
        endDisplacement = solidObject.qN1(nodeToMeasure, 1) - solidObject.qN(nodeToMeasure, 1);
        disp(['End displacement: ', num2str(endDisplacement)]);
        exactSolution = (1+solidObject.materialObject.nu)*p*3^2/(solidObject.materialObject.E*(9^2-3^2))*(9^2/evalPointX+(1-2*solidObject.materialObject.nu)*evalPointX);
        disp(['Ratio: ', num2str(endDisplacement/exactSolution)]);
    end
end
