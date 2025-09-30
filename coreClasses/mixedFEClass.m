classdef mixedFEClass < matlab.mixin.Copyable
    properties
        condensation = false;
        qR = [];
        qN = [];
        qN1 = [];
        globalEdof = [];
        shapeFunctionObject
        typeShapeFunction = 'sameOrder'; %'detailedOrder','other'
        typeShapeFunctionData = 0;
        continuousShapeFunctions = false;
        numberOfDofs = 0;
        numberOfFields = 0;
        dofsPerField = {};
        dofsPerFieldVector = [];
        dataFE = struct('ReDTilde',[],'KeDDTilde',[]);
    end
    methods
        function initializeMixedElements(obj,dofObject,continuumObject)
            obj.shapeFunctionObject = shapeFunctionClass();
            if strcmpi(obj.typeShapeFunction,'sameOrder')
                obj.shapeFunctionObject.order = obj.typeShapeFunctionData;
            end
            dimension = continuumObject.dimension;
            computeShapeFunctions = true;
            numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
            numberOfNodes = obj.shapeFunctionObject.numberOfNodes;
            switch continuumObject.elementDisplacementType
                case 'pianSumihara'
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.typeShapeFunctionData;
                    obj.numberOfDofs = obj.typeShapeFunctionData;
                    computeShapeFunctions = false;
                case {'incompatibleModesWilson','incompatibleModesTaylor'}
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = 4;
                    obj.numberOfDofs = 4;
                    computeShapeFunctions = false;
                case {'eas', 'easSC', 'incompatibleModes'}
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.typeShapeFunctionData;
                    obj.numberOfDofs = obj.typeShapeFunctionData;
                    computeShapeFunctions = false;
                    if isa(continuumObject,'solidVelocityClass')
                        obj.numberOfFields = 2;
                        obj.dofsPerField{1} = obj.typeShapeFunctionData;
                        obj.dofsPerField{2} = obj.typeShapeFunctionData;
                        obj.numberOfDofs = 2 * obj.typeShapeFunctionData;
                    end
                case {'incompressibleSimoTaylorPistor','incompressibleSimoTaylorPistorContinuous'}
                    if floor(obj.typeShapeFunctionData)~=ceil(obj.typeShapeFunctionData) || obj.typeShapeFunctionData < 0
                        error('invalid obj.orderShapeFunction')
                    end
                    obj.numberOfDofs = (obj.typeShapeFunctionData+1)^continuumObject.dimension*2;
                    obj.numberOfFields = 2;
                    obj.dofsPerField{1} = obj.numberOfDofs/2;
                    obj.dofsPerField{2} = obj.numberOfDofs/2;
                    %volumetric strain and pressure ansatz (lagrangian)
                    numberOfNodes = obj.numberOfDofs/2;
                case 'incompressibleSimoTaylorPistorBbar'
                    %check for legal ansatz
                    if obj.typeShapeFunctionData~=0 || continuumObject.shapeFunctionObject.order~=1
                        error('Simo-Taylor-Pister element as B-bar-method only implemented for Q1P0')
                    end
                    assert(obj.condensation, 'B-bar method can only be statically condensed');
                    % volumetric strain and pressure ansatz (lagrangian)
                    numberOfNodes = 1;
                case {'incompressibleHughes','incompressibleHerrmann'}
                    if floor(obj.typeShapeFunctionData)~=ceil(obj.typeShapeFunctionData) || obj.typeShapeFunctionData<0
                        error('invalid order of shapefunction for internal dof')
                    end
                    obj.numberOfDofs = (obj.typeShapeFunctionData+1)^continuumObject.dimension;
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.numberOfDofs;
                    %volumetric strain and pressure ansatz (lagrangian)
                    numberOfNodes = obj.numberOfDofs;
                case 'mixed'
                    obj.numberOfDofs = obj.typeShapeFunctionData;
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.numberOfDofs;
                case 'selectiveReducedIntegration'
                    assert(obj.condensation, 'selective reduced integration method can only be statically condensed');
                    numberOfGausspoints = (continuumObject.shapeFunctionObject.order)^dimension;
                    numberOfNodes = size(continuumObject.meshObject.edof,2);
                case {'mixedSC','mixedD_SC'}
                    switch continuumObject.elementDisplacementType
                        case {'mixedSC'}
                            if isa(continuumObject,'solidClass') || isa(continuumObject,'solidVelocityClass')
                                if (isempty(continuumObject.elementNameAdditionalSpecification) || strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeCGc') || strcmpi(continuumObject.elementNameAdditionalSpecification,'pHCGJLambda'))
                                    obj.numberOfFields = 6;
                                    obj.dofsPerFieldVector = [6;6;1;6;6;1];
                                elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'pHCGJ')
                                    obj.numberOfFields = 3;
                                    obj.dofsPerFieldVector = [6;6;1];
                                elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeGc') || strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeGJ') || strcmpi(continuumObject.elementNameAdditionalSpecification,'pHGJLambda')
                                    obj.numberOfFields = 4;
                                    obj.dofsPerFieldVector = [6;1;6;1];
                                elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'InvariantCGc') || strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeIIIJ')
                                    obj.numberOfFields = 6;
                                    obj.dofsPerFieldVector = [1;1;1;1;1;1];
                                end
                            elseif isa(continuumObject,'solidElectroClass') || isa(continuumObject,'solidElectroThermoClass')
                                obj.numberOfFields = 7;
                                obj.dofsPerFieldVector = [3;6;6;1;6;6;1];
                            else
                                error('class not implemented yet!');
                            end
                        case 'mixedD_SC'
                            obj.numberOfFields = 1;
                            obj.dofsPerFieldVector = 3;
                    end
                    if strcmpi(obj.typeShapeFunction,'sameOrder')
                        if isempty(numberOfNodes)
                            orderMixedShapeFunction = obj.typeShapeFunctionData;
                            if any(strcmpi(continuumObject.elementGeometryType, {'quadrilateral', 'hexahedral'}))
                                numberOfNodes = (orderMixedShapeFunction+1)^continuumObject.dimension;
                            elseif strcmpi(continuumObject.elementGeometryType, 'triangular')
                                if orderMixedShapeFunction == 0
                                    numberOfNodes = 1;
                                elseif orderMixedShapeFunction == 1
                                    numberOfNodes = 3;
                                elseif orderMixedShapeFunction == 2
                                    numberOfNodes = 6;
                                end
                            elseif strcmpi(continuumObject.elementGeometryType, 'tetrahedral')
                                if orderMixedShapeFunction == 0
                                    numberOfNodes = 1;
                                elseif orderMixedShapeFunction == 1
                                    numberOfNodes = 4;
                                elseif orderMixedShapeFunction == 2
                                    numberOfNodes = 10;
                                end
                            else
                                error('Element geometry type not implemented')
                            end
                        end
                        obj.dofsPerField = cell(obj.numberOfFields,1);
                        for ii = 1:obj.numberOfFields
                            obj.dofsPerField{ii,:} = [obj.dofsPerFieldVector(ii)*numberOfNodes obj.dofsPerFieldVector(ii)*numberOfNodes];
                        end
                        obj.numberOfDofs = sum(vertcat(obj.dofsPerField{:}),1);
                    elseif strcmpi(obj.typeShapeFunction,'detailedOrder')
                        obj.dofsPerField = cell(obj.numberOfFields,1);
                        numberOfNodes = zeros(obj.numberOfFields,2);
                        fieldsMixedShapeFunctions = fields(obj.typeShapeFunctionData);
                        for ii = 1:obj.numberOfFields
                            orderMixedShapeFunction = obj.typeShapeFunctionData.(fieldsMixedShapeFunctions{ii});
                            numberOfNodes(ii,:) = [(orderMixedShapeFunction(1)+1)^continuumObject.dimension (orderMixedShapeFunction(2)+1)^continuumObject.dimension];
                            obj.dofsPerField{ii,:} = [obj.dofsPerFieldVector(ii)*numberOfNodes(ii,1) obj.dofsPerFieldVector(ii)*numberOfNodes(ii,2)];
                        end
                        obj.numberOfDofs = sum(vertcat(obj.dofsPerField{:}),1);
                    elseif strcmpi(obj.typeShapeFunction,'other')
                        error('not implemented yet')
                    end
                case {'mixedPH'}
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = 1;
                    obj.numberOfDofs = 1; %extended strain variables
                    numberOfNodes = 1; % discontinuous ansatz constant elementwise
                    if isa(continuumObject,'beamClass')
                        obj.numberOfFields = 2;
                        if strcmpi(continuumObject.theory,'GeometricallyExact')
                            obj.dofsPerField{1} = 2;
                            obj.numberOfDofs = 3; % extended strain variables (elongation, shear deformation and curvature)
                        else
                            obj.dofsPerField{1} = 1;
                            obj.numberOfDofs = 2; % extended strain variables (shear deformation and curvature)
                        end
                        obj.dofsPerField{2} = 1;
                        numberOfNodes = 1; % discontinuous ansatz constant elementwise
                    elseif isa(continuumObject,'beamVelocityClass')
                        % TO CHECK: do sth else here?
                    end
                case {'mixedPHViscoPT'}
                    obj.numberOfFields = 2;
                    obj.dofsPerField{1} = 1;
                    obj.dofsPerField{2} = 1;
                    obj.numberOfDofs = 2; %extended strain variables and elastic strain
                    numberOfNodes = 1; % discontinuous ansatz constant elementwise
                case {'dispL2'}
                    obj.numberOfFields = 2;
                    obj.dofsPerFieldVector = [1; 1];

                    if isempty(numberOfNodes)
                        orderMixedShapeFunction = obj.typeShapeFunctionData;
                        if any(strcmpi(continuumObject.elementGeometryType, {'quadrilateral', 'hexahedral'}))
                            numberOfNodes = (orderMixedShapeFunction+1)^continuumObject.dimension;
                        elseif strcmpi(continuumObject.elementGeometryType, 'triangular')
                            if orderMixedShapeFunction == 0
                                numberOfNodes = 1;
                            elseif orderMixedShapeFunction == 1
                                numberOfNodes = 3;
                            elseif orderMixedShapeFunction == 2
                                numberOfNodes = 6;
                            end
                        elseif strcmpi(continuumObject.elementGeometryType, 'tetrahedral')
                            if orderMixedShapeFunction == 0
                                numberOfNodes = 1;
                            elseif orderMixedShapeFunction == 1
                                numberOfNodes = 4;
                            elseif orderMixedShapeFunction == 2
                                numberOfNodes = 10;
                            end
                        else
                            error('Element geometry type not implemented')
                        end
                    end
                    obj.dofsPerField = cell(obj.numberOfFields,1);
                    for ii = 1:obj.numberOfFields
                        obj.dofsPerField{ii,:} = [obj.dofsPerFieldVector(ii)*numberOfNodes obj.dofsPerFieldVector(ii)*numberOfNodes];
                    end
                    obj.numberOfDofs = sum(vertcat(obj.dofsPerField{:}),1);
                case {'mixedCGJL2'}
                    obj.numberOfFields = 5;
                    obj.dofsPerFieldVector = [1; 1; 6; 6; 1];

                    if isempty(numberOfNodes)
                        orderMixedShapeFunction = obj.typeShapeFunctionData;
                        if any(strcmpi(continuumObject.elementGeometryType, {'quadrilateral', 'hexahedral'}))
                            numberOfNodes = (orderMixedShapeFunction+1)^continuumObject.dimension;
                        elseif strcmpi(continuumObject.elementGeometryType, 'triangular')
                            if orderMixedShapeFunction == 0
                                numberOfNodes = 1;
                            elseif orderMixedShapeFunction == 1
                                numberOfNodes = 3;
                            elseif orderMixedShapeFunction == 2
                                numberOfNodes = 6;
                            end
                        elseif strcmpi(continuumObject.elementGeometryType, 'tetrahedral')
                            if orderMixedShapeFunction == 0
                                numberOfNodes = 1;
                            elseif orderMixedShapeFunction == 1
                                numberOfNodes = 4;
                            elseif orderMixedShapeFunction == 2
                                numberOfNodes = 10;
                            end
                        else
                            error('Element geometry type not implemented')
                        end
                    end
                    obj.dofsPerField = cell(obj.numberOfFields,1);
                    for ii = 1:obj.numberOfFields
                        obj.dofsPerField{ii,:} = [obj.dofsPerFieldVector(ii)*numberOfNodes obj.dofsPerFieldVector(ii)*numberOfNodes];
                    end
                    obj.numberOfDofs = sum(vertcat(obj.dofsPerField{:}),1);
                otherwise
                    error('elementDisplacementType not implemented')
            end

            % Compute mixed shape functions
            obj.shapeFunctionObject.numberOfNodes = numberOfNodes;
            if isempty(numberOfNodes)
                obj.shapeFunctionObject.numberOfNodes = continuumObject.shapeFunctionObject.numberOfNodes;
            end
            obj.shapeFunctionObject.numberOfGausspoints = numberOfGausspoints;
            if isempty(numberOfGausspoints)
                obj.shapeFunctionObject.numberOfGausspoints = continuumObject.shapeFunctionObject.numberOfGausspoints;
            end
            obj.shapeFunctionObject.computeMixedAnsatzFunctions(continuumObject, obj, computeShapeFunctions);


            numberOfElements = size(continuumObject.meshObject.edof,1);
            if obj.continuousShapeFunctions == false  % discontinous dofs
                if obj.condensation
                    % globalEdof (not needed if condensated)
                    obj.globalEdof = [];
                else
                    globalEdof = zeros(numberOfElements,obj.numberOfDofs(1))';
                    globalEdof(:) = dofObject.totalNumberOfDofs + (1:numberOfElements*obj.numberOfDofs(1));
                    obj.globalEdof = globalEdof';
                    dofObject.totalNumberOfDofs = dofObject.totalNumberOfDofs + numberOfElements*obj.numberOfDofs(1);
                end
            else % continuous dofs
                % globalEdof with connections
                if size(continuumObject.meshObject.edof,2) == 4 && obj.typeShapeFunctionData == 1
                    edofNodes = continuumObject.meshObject.edof;
                elseif size(continuumObject.meshObject.edof,2) == 9 && obj.typeShapeFunctionData == 1
                    edofNodes = continuumObject.meshObject.edof(:,1:4);
                else
                    error('combination of shapefunctions for continuous internal dofs is not implemented')
                end
                % squeeze the unused DOF
                uniqueDof = unique(edofNodes(:));
                for i = 1:size(edofNodes,1)
                    for j = 1:size(edofNodes,2)
                        % gets index of dof in uniqueDof lambda into edof --> local Dof without gaps
                        edofNodes(i,j) = find(edofNodes(i,j) == uniqueDof);
                    end
                end
                % adjust size for multiple Dof per node --> globalEdof
                dofPerNode = obj.numberOfDofs/size(edofNodes,2);
                numberOfNodes = size(edofNodes,2);
                if ceil(dofPerNode)~=floor(dofPerNode)
                    error('dofPerNode does not fit edof for continuous internal dof')
                end
                obj.globalEdof = zeros(numberOfElements,size(edofNodes,2)*dofPerNode);
                for j = 1:dofPerNode
                    obj.globalEdof(:,(j-1)*numberOfNodes+1:j*numberOfNodes) = dofObject.totalNumberOfDofs + edofNodes+(j-1)*size(uniqueDof,1);
                end
                dofObject.totalNumberOfDofs = dofObject.totalNumberOfDofs + size(uniqueDof,1)*dofPerNode;
            end
            % Set initial values for mixed fields
            if ~any(any(obj.qN))
                if numel(obj.numberOfDofs)==1
                    obj.qN = zeros(numberOfElements,obj.numberOfDofs);
                else
                    obj.qN = zeros(numberOfElements,obj.numberOfDofs(1));
                end
                if strcmp(continuumObject.elementDisplacementType, 'dispL2')
                    E = 1;
                    Eta = 1;
                    orderMixedShapeFunction = obj.typeShapeFunctionData;
                    numberOfNodes = (orderMixedShapeFunction+1)^continuumObject.dimension;
                    if numel(numberOfNodes)==1
                        obj.qN = [repmat(E,numberOfElements,numberOfNodes), repmat(Eta,numberOfElements,numberOfNodes)];
                    else
                        obj.qN = [repmat(E,numberOfElements,numberOfNodes(1)), repmat(Eta,numberOfElements,numberOfNodes(2))];
                    end
                elseif strcmp(continuumObject.elementDisplacementType, 'mixedCGJL2')
                    E = 1;
                    Eta = 1;
                    C   = eye(3);
                    CV  = matrixToVoigt(C,'strain').';
                    G = 1 / 2 * wedge(C,C);
                    GV = matrixToVoigt(G,'strain').';
                    J = sqrt(1 / 3 * innerProduct(G,C));
                    orderMixedShapeFunction = obj.typeShapeFunctionData;
                    numberOfNodes = (orderMixedShapeFunction+1)^continuumObject.dimension;
                    if numel(numberOfNodes)==1
                        obj.qN = [repmat(E,numberOfElements,numberOfNodes), repmat(Eta,numberOfElements,numberOfNodes), repmat(CV,numberOfElements,numberOfNodes), repmat(GV,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes)];
                    else
                        obj.qN = [repmat(E,numberOfElements,numberOfNodes(1)), repmat(Eta,numberOfElements,numberOfNodes(2)), repmat(CV,numberOfElements,numberOfNodes(3)), repmat(GV,numberOfElements,numberOfNodes(4)), repmat(J,numberOfElements,numberOfNodes(5))];
                    end
                elseif strcmp(continuumObject.elementDisplacementType, 'mixedSC')
                    if dimension == 3
                        switch continuumObject.materialObject.name
                            case {'ANN','MooneyRivlin', 'MooneyRivlinFullCoupled', 'MooneyRivlinFullCoupledMehnert','MooneyRivlinVol2','MooneyRivlinModified','MooneyRivlinModifiedVol2'}
                                C = eye(3);
                                Cv = matrixToVoigt(C, 'stress')';
                                G = 1/2*wedge(C, C);
                                Gv = matrixToVoigt(G, 'stress')';
                                c = 1/3*innerProduct(G, C);
                                %                                 numberOfNodes = size(obj.shapeFunctionObject.N, 2);
                                if isa(continuumObject,'solidClass') || isa(continuumObject,'solidVelocityClass')
                                    if isempty(continuumObject.elementNameAdditionalSpecification)
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_C = continuumObject.materialObject.a*eye(3);
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_c = -continuumObject.materialObject.d/(2*c)+continuumObject.materialObject.c/2*(1-1/sqrt(c));
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [dW_C, dW_G, dW_c, ~, ~, ~, ~, ~, ~] = ANNObject.computeDiffEnergyCGcANN(C,G,c);
                                        else
                                            error('material model not implemented')
                                        end

                                        lambdac = dW_c;
                                        lambdaG = dW_G + 1/3*lambdac*C;
                                        lambdaC = dW_C + wedge(lambdaG, C) + 1/3*lambdac*G;
                                        lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                        lambdaCv = matrixToVoigt(lambdaC, 'stress')';

                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(c,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdac,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes(1)), repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(c,numberOfElements,numberOfNodes(3)), repmat(lambdaCv,numberOfElements,numberOfNodes(4)), repmat(lambdaGv,numberOfElements,numberOfNodes(5)), repmat(lambdac,numberOfElements,numberOfNodes(6))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeCGc')
                                        C = eye(3);
                                        Cv = matrixToVoigt(C, 'strain')';
                                        G = 1/2*wedge(C, C);
                                        Gv = matrixToVoigt(G, 'strain')';
                                        c = 1/3*innerProduct(G, C);
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_C = continuumObject.materialObject.a*eye(3);
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_c = -continuumObject.materialObject.d/(2*c)+continuumObject.materialObject.c/2*(1-1/sqrt(c));
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [dW_C, dW_G, dW_c, ~, ~, ~, ~, ~, ~] = ANNObject.computeDiffEnergyCGcANN(C,G,c);
                                        else
                                            error('material model not implemented')
                                        end
                                        lambdac = dW_c;
                                        lambdaG = dW_G;
                                        lambdaC = dW_C;
                                        lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                        lambdaCv = matrixToVoigt(lambdaC, 'stress')';

                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(c,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdac,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes(1)), repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(c,numberOfElements,numberOfNodes(3)), repmat(lambdaCv,numberOfElements,numberOfNodes(4)), repmat(lambdaGv,numberOfElements,numberOfNodes(5)), repmat(lambdac,numberOfElements,numberOfNodes(6))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeGc')
                                        C = eye(3);
                                        G = 1/2*wedge(C, C);
                                        Gv = matrixToVoigt(G, 'strain')';
                                        c = 1/3*innerProduct(G, C);
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_II = continuumObject.materialObject.b;
                                            dW_III = -continuumObject.materialObject.d/(2*c)+continuumObject.materialObject.c/2*(1-1/sqrt(c));
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [~, ~, dW_III, ~, dW_II, ~, ~, ~, ~] = ANNObject.computeDiffEnergyCGcANN(C,G,c);
                                        else
                                            error('material model not implemented')
                                        end
                                        lambdac = dW_III;
                                        lambdaG = dW_II*eye(3);
                                        lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes), repmat(c,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdac,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes(1)), repmat(c,numberOfElements,numberOfNodes(2)), repmat(lambdaGv,numberOfElements,numberOfNodes(3)), repmat(lambdac,numberOfElements,numberOfNodes(4))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'pHCGJLambda')
                                        C = eye(3);
                                        G = eye(3);
                                        Cv = matrixToVoigt(C, 'stress')';
                                        Gv = matrixToVoigt(G, 'stress')';
                                        J = 1;
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_C = continuumObject.materialObject.a*eye(3);
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_J = continuumObject.materialObject.c*(J-1) - continuumObject.materialObject.d/J;
                                        elseif strcmpi(continuumObject.materialObject.name,'MooneyRivlinVol2')
                                            dW_C = continuumObject.materialObject.a*eye(3);
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_J = continuumObject.materialObject.c*(J-1) - continuumObject.materialObject.d;                                            
                                        elseif strcmpi(continuumObject.materialObject.name,'MooneyRivlinModified')
                                            alpha = continuumObject.materialObject.alpha;
                                            beta = continuumObject.materialObject.beta;
                                            gamma = continuumObject.materialObject.gamma;
                                            epsilon1 = continuumObject.materialObject.epsilon1;
                                            epsilon2 = continuumObject.materialObject.epsilon2;
                                            I = eye(3);
                                            dW_C = alpha * trace(C) * I;
                                            dW_G = beta * trace(G) * I;
                                            dW_J = - gamma*J^(-1) + 2*epsilon1*epsilon2*(J^(2*epsilon2-1) - J^(-2*epsilon2-1));
                                        end
                                        lambdaCv = matrixToVoigt(dW_C, 'stress')';
                                        lambdaGv = matrixToVoigt(dW_G, 'stress')';
                                        lambdaJ = dW_J;
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaJ,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes(1)), repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(J,numberOfElements,numberOfNodes(3)), repmat(lambdaCv,numberOfElements,numberOfNodes(4)), repmat(lambdaGv,numberOfElements,numberOfNodes(5)), repmat(lambdaJ,numberOfElements,numberOfNodes(6))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'pHGJLambda')
                                        G = eye(3);
                                        Gv = matrixToVoigt(G, 'stress')';
                                        J = 1;
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_J = continuumObject.materialObject.c*(J-1) - continuumObject.materialObject.d/J;
                                        elseif strcmpi(continuumObject.materialObject.name,'MooneyRivlinVol2')
                                            dW_G = continuumObject.materialObject.b*eye(3);
                                            dW_J = continuumObject.materialObject.c*(J-1) - continuumObject.materialObject.d;
                                        elseif strcmpi(continuumObject.materialObject.name,'MooneyRivlinModified')
                                            beta = continuumObject.materialObject.beta;
                                            gamma = continuumObject.materialObject.gamma;
                                            epsilon1 = continuumObject.materialObject.epsilon1;
                                            epsilon2 = continuumObject.materialObject.epsilon2;
                                            I = eye(3);
                                            dW_G = beta * trace(G) * I;
                                            dW_J = - gamma*J^(-1) + 2*epsilon1*epsilon2*(J^(2*epsilon2-1) - J^(-2*epsilon2-1));
                                        end
                                        lambdaGv = matrixToVoigt(dW_G, 'stress')';
                                        lambdaJ = dW_J;
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaJ,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(J,numberOfElements,numberOfNodes(3)), repmat(lambdaGv,numberOfElements,numberOfNodes(5)), repmat(lambdaJ,numberOfElements,numberOfNodes(6))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'pHCGJ')
                                        Cv = matrixToVoigt(eye(3), 'stress')';
                                        Gv = matrixToVoigt(eye(3), 'stress')';
                                        J = 1;
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Cv,numberOfElements,numberOfNodes(1)), repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(J,numberOfElements,numberOfNodes(3))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeGJ')
                                        C = eye(3);
                                        G = 1/2*wedge(C, C);
                                        Gv = matrixToVoigt(G, 'strain')';
                                        c = 1/3*innerProduct(G, C);
                                        J = sqrt(c);
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_II = continuumObject.materialObject.b;
                                            %                                             - d * log(J) + c / 2 * (J - 1)^2)
                                            dW_J = -continuumObject.materialObject.d/J + continuumObject.materialObject.c*(J-1);
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [~, dW_II, dW_J] = ANNObject.computeDiffEnergyANNJ(C,G,J);
                                        else
                                            error('material model not implemented')
                                        end
                                        lambdaJ = dW_J;
                                        lambdaG = dW_II*eye(3);
                                        lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaJ,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(Gv,numberOfElements,numberOfNodes(1)), repmat(J,numberOfElements,numberOfNodes(2)), repmat(lambdaGv,numberOfElements,numberOfNodes(3)), repmat(lambdaJ,numberOfElements,numberOfNodes(4))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'InvariantCGc')
                                        C = eye(3);
                                        I = trace(C);
                                        G = 1/2*wedge(C, C);
                                        II = trace(G);
                                        III = 1/3*innerProduct(G, C);
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_I = continuumObject.materialObject.a;
                                            dW_II = continuumObject.materialObject.b;
                                            dW_III = -continuumObject.materialObject.d / (2 * III) + continuumObject.materialObject.c / 2 * (1 - 1 / sqrt(III));
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [dW_I, dW_II, dW_III] = ANNObject.computeDiffEnergyInvariantANN(I,II,III);
                                        else
                                            error('material model not implemented')
                                        end
                                        lambdaIII = dW_III;
                                        lambdaII = dW_II;
                                        lambdaI = dW_I;
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(I,numberOfElements,numberOfNodes), repmat(II,numberOfElements,numberOfNodes), repmat(III,numberOfElements,numberOfNodes), repmat(lambdaI,numberOfElements,numberOfNodes), repmat(lambdaII,numberOfElements,numberOfNodes), repmat(lambdaIII,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(I,numberOfElements,numberOfNodes(1)), repmat(II,numberOfElements,numberOfNodes(2)), repmat(III,numberOfElements,numberOfNodes(3)), repmat(lambdaI,numberOfElements,numberOfNodes(4)), repmat(lambdaII,numberOfElements,numberOfNodes(5)), repmat(lambdaIII,numberOfElements,numberOfNodes(6))];
                                        end
                                    elseif strcmpi(continuumObject.elementNameAdditionalSpecification,'NonCascadeIIIJ')
                                        C = eye(3);
                                        I = trace(C);
                                        G = 1/2*wedge(C, C);
                                        II = trace(G);
                                        c = 1/3*innerProduct(G, C);
                                        J = sqrt(c);
                                        if strcmpi(continuumObject.materialObject.name,'MooneyRivlin')
                                            dW_I = continuumObject.materialObject.a;
                                            dW_II = continuumObject.materialObject.b;
                                            dW_J = -continuumObject.materialObject.d/J + continuumObject.materialObject.c*(J-1);
                                        elseif strcmpi(continuumObject.materialObject.name,'ANN')
                                            ANNObject = continuumObject.artificialNeuralNetworkObject;
                                            [dW_I, dW_II, dW_J] = ANNObject.computeDiffEnergyANNIIIJ(I,II,J);
                                        else
                                            error('material model not implemented')
                                        end
                                        lambdaI = dW_I;
                                        lambdaII = dW_II;
                                        lambdaJ = dW_J;
                                        if numel(numberOfNodes)==1
                                            obj.qN = [repmat(I,numberOfElements,numberOfNodes), repmat(II,numberOfElements,numberOfNodes), repmat(J,numberOfElements,numberOfNodes), repmat(lambdaI,numberOfElements,numberOfNodes), repmat(lambdaII,numberOfElements,numberOfNodes), repmat(lambdaJ,numberOfElements,numberOfNodes)];
                                        else
                                            obj.qN = [repmat(I,numberOfElements,numberOfNodes(1)), repmat(II,numberOfElements,numberOfNodes(2)), repmat(J,numberOfElements,numberOfNodes(3)), repmat(lambdaI,numberOfElements,numberOfNodes(4)), repmat(lambdaII,numberOfElements,numberOfNodes(5)), repmat(lambdaJ,numberOfElements,numberOfNodes(6))];
                                        end
                                    end
                                elseif isa(continuumObject,'solidElectroClass')
                                    D = [0, 0, 0];
                                    dW_C = continuumObject.materialObject.a*eye(3);
                                    dW_G = continuumObject.materialObject.b*eye(3);
                                    dW_c = -continuumObject.materialObject.d/(2*c)+continuumObject.materialObject.c/2*(1-1/sqrt(c));

                                    lambdac = dW_c;
                                    lambdaG = dW_G + 1/3*lambdac*C;
                                    lambdaC = dW_C + wedge(lambdaG, C) + 1/3*lambdac*G;
                                    lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                    lambdaCv = matrixToVoigt(lambdaC, 'stress')';
                                    if numel(numberOfNodes)==1
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes), repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(c,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdac,numberOfElements,numberOfNodes)];
                                    else
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes(1)), repmat(Cv,numberOfElements,numberOfNodes(2)), repmat(Gv,numberOfElements,numberOfNodes(3)), repmat(c,numberOfElements,numberOfNodes(4)), repmat(lambdaCv,numberOfElements,numberOfNodes(5)), repmat(lambdaGv,numberOfElements,numberOfNodes(6)), repmat(lambdac,numberOfElements,numberOfNodes(7))];
                                    end
                                elseif isa(continuumObject,'solidElectroThermoClass')
                                    D = [0, 0, 0];
                                    % Note: works only for a constant temperature distribution in the initial configuration
                                    theta = sum(continuumObject.qN(:, continuumObject.dimension + 2)) / size(continuumObject.qN, 1);
                                    thetaR = continuumObject.materialObject.thetaR;
                                    beta = continuumObject.materialObject.beta;
                                    c2 = continuumObject.materialObject.c2;
                                    if strcmpi(continuumObject.materialObject.name, 'MooneyRivlin')
                                        dW_C = continuumObject.materialObject.a*eye(3);
                                        dW_G = continuumObject.materialObject.b*eye(3);
                                        dW_c = -continuumObject.materialObject.d1/(2*c)+continuumObject.materialObject.c1/2*(1-1/sqrt(c)) - 3*beta*c2*(theta - thetaR);
                                    elseif strcmpi(continuumObject.materialObject.name, 'MooneyRivlinFullCoupled')
                                        dW_C = theta/thetaR *continuumObject.materialObject.a*eye(3);
                                        dW_G = theta/thetaR * continuumObject.materialObject.b*eye(3);
                                        dW_c = theta/thetaR *(-continuumObject.materialObject.d1/(2*c)+continuumObject.materialObject.c1/2*(1-1/sqrt(c))) - 3*beta*c2*(theta - thetaR);
                                    end
                                    lambdac = dW_c;
                                    lambdaG = dW_G + 1/3*lambdac*C;
                                    lambdaC = dW_C + wedge(lambdaG, C) + 1/3*lambdac*G;
                                    lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                    lambdaCv = matrixToVoigt(lambdaC, 'stress')';
                                    if numel(numberOfNodes)==1
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes), repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(c,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdac,numberOfElements,numberOfNodes)];
                                    else
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes(1)), repmat(Cv,numberOfElements,numberOfNodes(2)), repmat(Gv,numberOfElements,numberOfNodes(3)), repmat(c,numberOfElements,numberOfNodes(4)), repmat(lambdaCv,numberOfElements,numberOfNodes(5)), repmat(lambdaGv,numberOfElements,numberOfNodes(6)), repmat(lambdac,numberOfElements,numberOfNodes(7))];
                                    end
                                else
                                    error('mixedSC class not implemented yet!');
                                end
                            otherwise
                                error('materialObject.name not implemented')
                        end
                    end
                elseif strcmp(continuumObject.elementDisplacementType, 'mixedD_SC')
                    if dimension == 3
                        switch continuumObject.materialObject.name
                            case {'MooneyRivlin'}
                                % not neccessary
                            otherwise
                                error('materialObject.name not implemented')
                        end
                    end
                elseif strcmp(continuumObject.elementDisplacementType, 'mixedPH')
                    if strcmp(continuumObject.elementNameAdditionalSpecification,'C')
                        if isempty(obj.qR) %no initial values for mixed quantity specified
                            obj.qN = ones(numberOfElements,obj.numberOfDofs);
                        else
                            obj.qN = obj.qR;
                        end
                    elseif strcmp(continuumObject.elementNameAdditionalSpecification,'E')
                        obj.qN = zeros(numberOfElements,obj.numberOfDofs);
                    end
                elseif strcmp(continuumObject.elementDisplacementType, 'mixedPHViscoPT')

                    if strcmp(continuumObject.elementNameAdditionalSpecification,'C')

                        if isempty(obj.qR) %no reference values for mixed quantity specified
                            obj.qN = ones(numberOfElements,obj.numberOfDofs);
                        else
                            obj.qN = obj.qR;
                        end
                    end
                end
            end

            obj.qN1 = obj.qN;
        end
        function updateGlobalField(obj,continuumObject,dofObject)
            % for preinitialization of internalDofs/mixedFields without condensation for dofObject.qN array
            if contains(continuumObject.elementDisplacementType, 'mixed')
                % if ~contains(continuumObject.elementDisplacementType, 'displacement')
                % does not work for all other mixed elements i.e. for Q1Theta1P1discont
                if ~obj.condensation
                    numberOfElements = size(continuumObject.meshObject.globalFullEdof, 1);
                    numberOfDisplacementDofs = dofObject.totalNumberOfDofs - numberOfElements*obj.numberOfDofs(1);
                    internalDofs = obj.qN';
                    dofObject.qN = dofObject.qN + [zeros(numberOfDisplacementDofs,1); internalDofs(:)];
                end
            end
        end
        function [array, dataFEMixed] = mixedFEFunction(obj,setupObject,RD,RA,KDD,KDA,KAD,KAA,varargin)
            if ~isempty(varargin)
                % pivotization/permutations
                P = varargin{1};
                RD = P*RD;
                KDD = P*KDD*P';
                if obj.numberOfFields >= 1
                    KDA = P*KDA;
                    KAD = KAD*P';
                end
            end
            if obj.condensation && obj.numberOfFields >= 1
                ReDTilde = solverLinearSystem(setupObject, KAA, RA);
                array.Re  = RD - KDA*ReDTilde;
                KeDDTilde = solverLinearSystem(setupObject, KAA, KAD);
                array.Ke  = KDD - KDA*KeDDTilde;
                dataFEMixed.ReDTilde = ReDTilde;
                dataFEMixed.KeDDTilde = KeDDTilde;
            else
                array.Re = [RD; RA];
                array.Ke = [KDD, KDA; KAD, KAA];
                dataFEMixed.ReDTilde = [];
                dataFEMixed.KeDDTilde = [];
            end
        end
        function updateMixedElements(obj,continuumObject,dofObject,setupObject,delta,fieldName,field)
            if obj.condensation == false
                if size(obj.globalEdof,1) == 1
                    obj.(fieldName) = field(obj.globalEdof)';
                else
                    obj.(fieldName) = field(obj.globalEdof);
                end
            else
                if any(contains(continuumObject.elementDisplacementType,{'eas', 'mixedSC', 'pianSumihara', 'mixedD_SC'}))
                    deltaQ = zeros(size(obj.(fieldName)));
                    for e = 1:size(continuumObject.meshObject.globalFullEdof,1)
                        if setupObject.newton.step == 1
                            delta = dofObject.qN1-dofObject.qN;
                        end
                        deltaQ(e,:) = -(obj.dataFE(e).ReDTilde + obj.dataFE(e).KeDDTilde*delta(continuumObject.meshObject.globalFullEdof(e,:)))';
                    end
                    obj.(fieldName) = obj.(fieldName) + deltaQ;
                end
            end
        end
    end
end
