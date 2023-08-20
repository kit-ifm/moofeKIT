classdef shapeFunctionClass < matlab.mixin.Copyable
    properties
        numberOfGausspoints
        order % can be removed?!
        numberOfNodes
    end
    properties (SetAccess = private)
        dimension
        gaussWeight % gauss weights
        gaussPoint % gauss points
        N_k_I % Lagrange shape functions
        dN_xi_k_I % First derivative of Lagrange shape functions
        d2N_xi_xi_k_I % Second derivative of Lagrange shape functions
        N0_k_I % Lagrange shape functions evaluated at the element's center
        dN0_xi_I % First derivative of Lagrange shape functions evaluated at the element's center

        M % Other ansatz functions (e.g. EAS, Pian-Sumihara, ...)
        dM % First derivative of other ansatz functions
    end
    methods
        function computeShapeFunction(obj, dimension, numberOfNodesBasedOnMeshObject, elementGeometryType)
            % set dimension
            obj.dimension = dimension;

            % compute gauss points and weights
            [obj.gaussPoint, obj.gaussWeight] = gaussPointsAndWeights(dimension, obj.numberOfGausspoints, elementGeometryType);

            % determine numberOfNodesArray
            numberOfNodesArray = obj.numberOfNodes;
            if isempty(numberOfNodesArray) && length(numberOfNodesBasedOnMeshObject) == 1
                numberOfNodesArray = numberOfNodesBasedOnMeshObject;
            end
            numberOfShapeFunctions = length(numberOfNodesArray);

            % compute shape functions
            NCell = cell(numberOfShapeFunctions, 1);
            dNCell = cell(numberOfShapeFunctions, 1);
            d2NCell = cell(numberOfShapeFunctions, 1);
            N0Cell = cell(numberOfShapeFunctions, 1);
            dN0Cell = cell(numberOfShapeFunctions, 1);
            for ii = 1:numberOfShapeFunctions
                % Lagrange shape functions evaluated at gauss points
                [NCell{ii}, dNCell{ii}, d2NCell{ii}] = computeLagrangeShapeFunction(dimension, numberOfNodesArray(ii), obj.numberOfGausspoints, obj.gaussPoint);

                % Lagrange shape functions evaluated at the elements center
                [gaussPoints, ~] = gaussPointsAndWeights(dimension, 1, elementGeometryType);
                [N0Cell{ii}, dN0_xi_k_I, ~] = computeLagrangeShapeFunction(dimension, numberOfNodesArray(ii), 1, gaussPoints);
                dN0Cell{ii} = reshape(dN0_xi_k_I(:, 1, :), [size(dN0_xi_k_I, 1), size(dN0_xi_k_I, 3)]);
            end

            % set values
            if numberOfShapeFunctions == 1
                obj.N_k_I = NCell{1};
                obj.dN_xi_k_I = dNCell{1};
                obj.d2N_xi_xi_k_I = d2NCell{1};
                obj.N0_k_I = N0Cell{1};
                obj.dN0_xi_I = dN0Cell{1};
            else
                obj.N_k_I = NCell;
                obj.dN_xi_k_I = dNCell;
                obj.d2N_xi_xi_k_I = d2NCell;
                obj.N0_k_I = N0Cell;
                obj.dN0_xi_I = dN0Cell;
            end
            obj.numberOfNodes = numberOfNodesArray;
        end
        function computeShapeFunctionCellArray(obj, dimension, numberOfNodes, elementGeometryType, typeShapeFunctionData)
            obj.dimension = dimension;
            obj.numberOfNodes = numberOfNodes;
            for ii = 1:numel(fieldnames(typeShapeFunctionData))
                for jj = 1:2
                    [obj.gaussPoint, obj.gaussWeight] = gaussPointsAndWeights(dimension, obj.numberOfGausspoints, elementGeometryType);
                    [N_k_I, dN_xi_k_I, d2N_xi_k_I] = computeLagrangeShapeFunction(dimension, numberOfNodes(ii, jj), obj.numberOfGausspoints, obj.gaussPoint);
                    obj.N_k_I{ii, jj} = N_k_I;
                    obj.dN_xi_k_I{ii, jj} = dN_xi_k_I;
                    obj.dNxi2{ii, jj} = d2N_xi_k_I;
                end
            end
        end

        function computeMixedAnsatzFunctions(obj, continuumObject, mixedFEObject, computeShapeFunctions)
            % compute Lagrange shape functions and gauss points
            if strcmpi(mixedFEObject.typeShapeFunction,'detailedOrder')
                assert(mixedFEObject.numberOfFields == numel(fieldnames(mixedFEObject.typeShapeFunctionData)), 'number of fields did not match');
                obj.computeShapeFunctionCellArray(continuumObject.dimension, obj.numberOfNodes, continuumObject.elementGeometryType, obj.typeShapeFunctionData);
            else
                obj.computeShapeFunction(continuumObject.dimension, obj.numberOfNodes, continuumObject.elementGeometryType);
            end

            % compute mixed ansatz functions
            switch continuumObject.elementDisplacementType
                case 'pianSumihara'
                    [obj.M, ~] = computePianSumiharaAnsatzFunctions(continuumObject, obj.dimension, mixedFEObject.typeShapeFunctionData, obj.numberOfGausspoints, obj.gaussPoint);
                case {'incompatibleModesWilson','incompatibleModesTaylor'}
                    %bubble modes
                    M_k_I = zeros(obj.numberOfGausspoints, obj.dimension);
                    dM_xi_k_I = zeros(obj.dimension, obj.numberOfGausspoints, obj.dimension);
                    M_k_I(:,1) = 1-obj.gaussPoint(1,:).^2;
                    M_k_I(:,2) = 1-obj.gaussPoint(2,:).^2;
                    %derivatives w.r.t. x
                    dM_xi_k_I(1,:,1) = -2*obj.gaussPoint(1,:);
                    dM_xi_k_I(1,:,2) = 0;
                    %derivatives w.r.t. y
                    dM_xi_k_I(2,:,1) = 0;
                    dM_xi_k_I(2,:,2) = -2*obj.gaussPoint(2,:);
                    obj.M = M_k_I;
                    obj.dM = dM_xi_k_I;
                case {'eas', 'easSC', 'easPetrovGalerkin'}
                    [obj.M, ~] = computeEASAnsatzFunctions(continuumObject, obj.dimension, mixedFEObject.typeShapeFunctionData, obj.numberOfGausspoints, obj.gaussPoint);
            end
        end
    end

    % set and get methods
    methods
        function set.dimension(obj, input)
            checkIntegerMatrix(input);
            assert((input <= 3) & (input >= 1));
            obj.dimension = input;
        end
    end
end