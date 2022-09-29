classdef mixedFEClass < matlab.mixin.Copyable
    properties
        condensation = false;
        qN = [];
        qN1 = [];
        globalEdof = [];
        shapeFunctionObject
        typeShapeFunction = 'sameOrder'; %'detailedOrder','other'
        typeShapeFunctionData = 0;
        continuousShapeFunctions = false;
        numberOfNodes
        numberOfDofs = 0;
        numberOfFields = 0;
        dofsPerField = {};
        dofsPerFieldVector = [];
        dataFE = struct('ReDTilde',[],'KeDDTilde',[]);
    end
    methods
        function initializeMixedElements(obj,dofObject,continuumObject)
            obj.shapeFunctionObject = shapeFunctionClass();
            shapeFunctionObject = obj.shapeFunctionObject;
            if strcmpi(obj.typeShapeFunction,'sameOrder')
                shapeFunctionObject.order = obj.typeShapeFunctionData;
            end
            dimension = continuumObject.dimension;
            numberOfGausspoints = continuumObject.shapeFunctionObject.numberOfGausspoints;
            gaussPoint = continuumObject.shapeFunctionObject.gaussPoint';
            computeShapeFunctions = true;
            centerEvaluationDisplacementShapeFunction = false;
            switch continuumObject.elementDisplacementType
                case 'pianSumihara'
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.typeShapeFunctionData;
                    obj.numberOfDofs = obj.typeShapeFunctionData;
                    centerEvaluationDisplacementShapeFunction = true;
                    [N_k_I, dN_xi_k_I] = computePianSumiharaAnsatzFunctions(dimension, obj.typeShapeFunctionData, numberOfGausspoints, gaussPoint);
                    computeShapeFunctions = false;
                case {'incompatibleModesWilson','incompatibleModesTaylor'}
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = 4;
                    obj.numberOfDofs = 4;
                    centerEvaluationDisplacementShapeFunction = true;
                    %bubble modes
                    N_k_I = zeros(numberOfGausspoints,dimension);
                    dN_xi_k_I = zeros(dimension,numberOfGausspoints,dimension);
                    N_k_I(:,1) = 1-gaussPoint(:,1).^2;
                    N_k_I(:,2) = 1-gaussPoint(:,2).^2;
                    %derivatives w.r.t. x
                    dN_xi_k_I(1,:,1) = -2*gaussPoint(:,1);
                    dN_xi_k_I(1,:,2) = 0;
                    %derivatives w.r.t. y
                    dN_xi_k_I(2,:,1) = 0;
                    dN_xi_k_I(2,:,2) = -2*gaussPoint(:,2);
                    computeShapeFunctions = false;
                case {'eas', 'easSC', 'easPetrovGalerkin'}
                    obj.numberOfFields = 1;
                    obj.dofsPerField{1} = obj.typeShapeFunctionData;
                    obj.numberOfDofs = obj.typeShapeFunctionData;
                    [N_k_I, dN_xi_k_I] = computeEASAnsatzFunctions(dimension, obj.typeShapeFunctionData, numberOfGausspoints, gaussPoint);
                    centerEvaluationDisplacementShapeFunction = true;
                    computeShapeFunctions = false;
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
                    if obj.condensation==false
                        error('B-bar method can only be statically condensed')
                    end
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
                case 'selectiveReducedIntegration'
                    if obj.condensation==false
                        error('selective reduced integration method can only be statically condensed')
                    end
                    numberOfGausspoints = (continuumObject.shapeFunctionObject.order)^dimension;
                    numberOfNodes = size(continuumObject.meshObject.edof,2);
                case {'mixedSC','mixedD_SC'}
                    switch continuumObject.elementDisplacementType
                        case 'mixedSC'
                            if isa(continuumObject,'solidClass')
                                obj.numberOfFields = 6;
                                obj.dofsPerFieldVector = [6;6;1;6;6;1];
                            elseif isa(continuumObject,'solidElectroClass') || isa(continuumObject,'solidElectroThermoClass')
                                obj.numberOfFields = 7;
                                obj.dofsPerFieldVector = [3;6;6;1;6;6;1];
                            else
                                error('class not implemented yet!');
                            end
                        case 'mixedD_SC'
                            obj.numberOfFields = 1;
                            obj.dofsPerFieldVector = [3];
                    end
                    if strcmpi(obj.typeShapeFunction,'sameOrder')
                        if ~isempty(obj.numberOfNodes)
                            numberOfNodes = obj.numberOfNodes;
                        else
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
                otherwise
                    error('elementDisplacementType not implemented')
            end
            if computeShapeFunctions
                shapeFunctionObject.numberOfGausspoints = numberOfGausspoints;
                if strcmpi(obj.typeShapeFunction,'detailedOrder')
                    assert(obj.numberOfFields == numel(fieldnames(obj.typeShapeFunctionData)),'number of fields did not match');
                    shapeFunctionObject.computeShapeFunctionCellArray(dimension,numberOfNodes,continuumObject.elementGeometryType,obj.typeShapeFunctionData);
                else
                    shapeFunctionObject.computeShapeFunction(dimension,numberOfNodes,continuumObject.elementGeometryType);
                end
            else
                shapeFunctionObject.dN_xi_k_I = dN_xi_k_I;
                shapeFunctionObject.N_k_I = N_k_I;
            end
            if centerEvaluationDisplacementShapeFunction
                numberOfGausspoints = 1;
                [gaussPoints, ~] = gaussPointsAndWeights(dimension, numberOfGausspoints, continuumObject.elementGeometryType);
                [~,dNxi0,~] = computeLagrangeShapeFunction(dimension,continuumObject.shapeFunctionObject.numberOfNodes,numberOfGausspoints,gaussPoints);
                continuumObject.shapeFunctionObject.dN0_xi_k_I = dNxi0;
            end
            numberOfElements = size(continuumObject.meshObject.edof,1);
            if obj.continuousShapeFunctions == false  % discontinous dofs
                if obj.condensation == true 
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
                if strcmp(continuumObject.elementDisplacementType, 'mixedSC')
                    if dimension == 3
                        switch continuumObject.materialObject.name
                            case {'MooneyRivlin', 'MooneyRivlinFullCoupled', 'MooneyRivlinFullCoupledMehnert'}
                                C = eye(3);
                                Cv = matrixToVoigt(C, 'stress')';
                                G = 1/2*wedge(C, C);
                                Gv = matrixToVoigt(G, 'stress')';
                                I3 = 1/3*innerProduct(G, C);
%                                 numberOfNodes = size(obj.shapeFunctionObject.N, 2);
                                if isa(continuumObject,'solidClass')
                                    Sigma_C = continuumObject.materialObject.a*eye(3);
                                    Sigma_G = continuumObject.materialObject.b*eye(3);
                                    Sigma_I = -continuumObject.materialObject.d/(2*I3)+continuumObject.materialObject.c/2*(1-1/sqrt(I3));

                                    lambdaI3 = Sigma_I;
                                    lambdaG = Sigma_G + 1/3*lambdaI3*C;
                                    lambdaC = Sigma_C + wedge(lambdaG, C) + 1/3*lambdaI3*G;
                                    lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                    lambdaCv = matrixToVoigt(lambdaC, 'stress')';

                                    if numel(numberOfNodes)==1
                                        obj.qN = [repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(I3,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaI3,numberOfElements,numberOfNodes)];
                                    else
                                        obj.qN = [repmat(Cv,numberOfElements,numberOfNodes(1)), repmat(Gv,numberOfElements,numberOfNodes(2)), repmat(I3,numberOfElements,numberOfNodes(3)), repmat(lambdaCv,numberOfElements,numberOfNodes(4)), repmat(lambdaGv,numberOfElements,numberOfNodes(5)), repmat(lambdaI3,numberOfElements,numberOfNodes(6))];
                                    end                                    
                                elseif isa(continuumObject,'solidElectroClass')
                                    D = [0, 0, 0];
                                    Sigma_C = continuumObject.materialObject.a*eye(3);
                                    Sigma_G = continuumObject.materialObject.b*eye(3);
                                    Sigma_I = -continuumObject.materialObject.d/(2*I3)+continuumObject.materialObject.c/2*(1-1/sqrt(I3));

                                    lambdaI3 = Sigma_I;
                                    lambdaG = Sigma_G + 1/3*lambdaI3*C;
                                    lambdaC = Sigma_C + wedge(lambdaG, C) + 1/3*lambdaI3*G;
                                    lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                    lambdaCv = matrixToVoigt(lambdaC, 'stress')';
                                    if numel(numberOfNodes)==1
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes), repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(I3,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaI3,numberOfElements,numberOfNodes)];
                                    else
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes(1)), repmat(Cv,numberOfElements,numberOfNodes(2)), repmat(Gv,numberOfElements,numberOfNodes(3)), repmat(I3,numberOfElements,numberOfNodes(4)), repmat(lambdaCv,numberOfElements,numberOfNodes(5)), repmat(lambdaGv,numberOfElements,numberOfNodes(6)), repmat(lambdaI3,numberOfElements,numberOfNodes(7))];
                                    end
                                elseif isa(continuumObject,'solidElectroThermoClass')
                                    D = [0, 0, 0];
                                    % Note: works only for a constant temperature distribution in the initial configuration
                                    theta = sum(continuumObject.qN(:, continuumObject.dimension + 2)) / size(continuumObject.qN, 1);
                                    thetaR = continuumObject.materialObject.thetaR;
                                    beta = continuumObject.materialObject.beta;
                                    c2 = continuumObject.materialObject.c2;
                                    if strcmpi(continuumObject.materialObject.name, 'MooneyRivlin')
                                        Sigma_C = continuumObject.materialObject.a*eye(3);
                                        Sigma_G = continuumObject.materialObject.b*eye(3);
                                        Sigma_I = -continuumObject.materialObject.d1/(2*I3)+continuumObject.materialObject.c1/2*(1-1/sqrt(I3)) - 3*beta*c2*(theta - thetaR);
                                    elseif strcmpi(continuumObject.materialObject.name, 'MooneyRivlinFullCoupled')
                                        Sigma_C = theta/thetaR *continuumObject.materialObject.a*eye(3);
                                        Sigma_G = theta/thetaR * continuumObject.materialObject.b*eye(3);
                                        Sigma_I = theta/thetaR *(-continuumObject.materialObject.d1/(2*I3)+continuumObject.materialObject.c1/2*(1-1/sqrt(I3))) - 3*beta*c2*(theta - thetaR);
                                    end
                                    lambdaI3 = Sigma_I;
                                    lambdaG = Sigma_G + 1/3*lambdaI3*C;
                                    lambdaC = Sigma_C + wedge(lambdaG, C) + 1/3*lambdaI3*G;
                                    lambdaGv = matrixToVoigt(lambdaG, 'stress')';
                                    lambdaCv = matrixToVoigt(lambdaC, 'stress')';
                                    if numel(numberOfNodes)==1
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes), repmat(Cv,numberOfElements,numberOfNodes), repmat(Gv,numberOfElements,numberOfNodes), repmat(I3,numberOfElements,numberOfNodes), repmat(lambdaCv,numberOfElements,numberOfNodes), repmat(lambdaGv,numberOfElements,numberOfNodes), repmat(lambdaI3,numberOfElements,numberOfNodes)];
                                    else
                                        obj.qN = [repmat(D,numberOfElements,numberOfNodes(1)), repmat(Cv,numberOfElements,numberOfNodes(2)), repmat(Gv,numberOfElements,numberOfNodes(3)), repmat(I3,numberOfElements,numberOfNodes(4)), repmat(lambdaCv,numberOfElements,numberOfNodes(5)), repmat(lambdaGv,numberOfElements,numberOfNodes(6)), repmat(lambdaI3,numberOfElements,numberOfNodes(7))];
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
                end
            end
            obj.qN1 = obj.qN;
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
