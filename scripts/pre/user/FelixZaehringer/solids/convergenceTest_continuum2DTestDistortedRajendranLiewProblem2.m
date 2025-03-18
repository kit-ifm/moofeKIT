% run('../../../startUpMoofeKIT.m');

%% Start parameters
elementsToTest = {'disp', 'dispPetrovGalerkin', 'dispQ9'}; % {'disp', 'BD', 'eas', 'easPetrovGalerkin', 'easPetrovGalerkinBD', 'dispPetrovGalerkinXie'};
numberOfRefinements = 8;
maxDistortion = 25;

%% Script
resultArray = zeros(numberOfRefinements+1, size(elementsToTest, 2));
hArray = zeros(numberOfRefinements+1, 1);
numberOfDofArray = zeros(numberOfRefinements+1, 1);
numberOfElementsArray = zeros(numberOfRefinements+1, 1);
for ii = 1:size(elementsToTest, 2)
    for jj = 1:numberOfRefinements+1
%         s = (jj - 1)*maxDistortion/numberStepsMeshDistortion;

        lengthX = 100;
        lengthY = 5;
        h = lengthX/jj;
        hArray(jj) = h;

        disp('=======================================');
        disp(['Current element: ', elementsToTest{ii}]);
        disp(['Current h: ', num2str(h)]);
        disp('=======================================');

        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.saveObject.fileName = 'continuum2DTestDistorted';
        setupObject.saveObject.saveData = false;
        setupObject.totalTimeSteps = 1;
        setupObject.totalTime = 1;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = 11;
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

        [solidObject.meshObject.nodes,solidObject.meshObject.edof, bounEdof] = meshRectangle(lengthX, lengthY, jj*2, jj*2, order, serendipity);
%         solidObject.meshObject.nodes = distortedMeshGenerator(solidObject.meshObject.nodes, lengthX, lengthY);
        solidObject.meshObject.nodes = adaptMeshToCornerNodes(solidObject.meshObject.nodes, solidObject.meshObject.edof, 2);
        solidObject.meshObject.nodes(:,1) = solidObject.meshObject.nodes(:,1) + lengthX/2;
        solidObject.meshObject.nodes(:,2) = solidObject.meshObject.nodes(:,2) + lengthY/2;

        numberOfDofArray(jj) = length(unique(solidObject.meshObject.edof));
        numberOfElementsArray (jj) = size(solidObject.meshObject.edof, 1);

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
        
        % Neumann boundary
        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) [0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX2;

        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) -[0; 120*XYZ(2)/lengthX-120*XYZ(2)^2/(lengthX*lengthY)]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX1;

        neumannObject = neumannClass(dofObject);
        neumannObject.loadGeometry = 'line';
        neumannObject.masterObject = solidObject;
        neumannObject.loadVectorFunction = @(XYZ) [240*XYZ(2)/lengthY-120; 0]; % Komponenten des Kraftvektors anpassen
        neumannObject.meshObject.edof = bounEdof.SX1;
        
        %% solver
        dofObject = runNewton(setupObject,dofObject);
        evalPointX = lengthX;
        evalPointY = 0;
        nodeToMeasure = find(solidObject.meshObject.nodes(:,1) == evalPointX & solidObject.meshObject.nodes(:,2) == evalPointY);
        endDisplacement = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
        exactSolution = (-40/(lengthX*lengthY)*evalPointX^3-36/(lengthX*lengthY)*evalPointX*evalPointY^2+120/(lengthY)*evalPointX^2+36/lengthX*evalPointX*evalPointY+36/lengthY*evalPointY^2+46*lengthY/lengthX*evalPointX-36*evalPointY)/solidObject.materialObject.E;
        disp(['End displacement: ', num2str(endDisplacement)]);
        error = abs(exactSolution-endDisplacement)/abs(exactSolution);
        disp(['Error: ', num2str(error)]);
        resultArray(jj, ii) = error;
        
    end
end

evalStart = 1;
xVals = numberOfDofArray;
figure;
for ii=1:size(elementsToTest, 2)
    fit=polyfit(log(xVals(evalStart:end)),log(resultArray(evalStart:end, ii)),1);
    x1=min(xVals(evalStart:end)):0.1:max(xVals(evalStart:end));
    y1=exp(fit(1)*log(x1)+fit(2)); 
%     loglog(xVals(evalStart:end), resultArray(evalStart:end, ii), 'o', x1, y1, '-');
    loglog(xVals(evalStart:end), resultArray(evalStart:end, ii), 'o');
    hold on;
    disp(['Steigung: ', num2str(fit(1))])
end

xlabel('s');
ylabel('w_{max}');
legend(elementsToTest);
