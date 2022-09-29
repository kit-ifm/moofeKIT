classdef shapeFunctionClass < matlab.mixin.Copyable
    properties
        numberOfGausspoints
        order
    end
    properties
        gaussWeight
        gaussPoint
        numberOfNodes
        dimension
        N_k_I
        dN_xi_k_I
        d2N_xi_xi_k_I
        dN0_xi_k_I
    end
    methods
        function computeShapeFunction(obj, dimension, numberOfNodes, elementGeometryType)
            obj.dimension = dimension;
            obj.numberOfNodes = numberOfNodes;
            [obj.gaussPoint, obj.gaussWeight] = gaussPointsAndWeights(dimension, obj.numberOfGausspoints, elementGeometryType);
            [obj.N_k_I, obj.dN_xi_k_I, obj.d2N_xi_xi_k_I] = computeLagrangeShapeFunction(dimension, numberOfNodes, obj.numberOfGausspoints, obj.gaussPoint);
        end
        function computeShapeFunctionCellArray(obj,dimension, numberOfNodes, elementGeometryType, typeShapeFunction)
            obj.dimension = dimension;
            obj.numberOfNodes = numberOfNodes;
            for ii = 1:numel(fieldnames(typeShapeFunction))
                for jj = 1:2
                    [obj.gaussPoint, obj.gaussWeight] = gaussPointsAndWeights(dimension, obj.numberOfGausspoints, elementGeometryType);
                    [N_k_I, dN_xi_k_I, d2N_xi_k_I] = computeLagrangeShapeFunction(dimension,numberOfNodes(ii,jj),obj.numberOfGausspoints, obj.gaussPoint);
                    obj.N_k_I{ii,jj} = N_k_I;
                    obj.dN_xi_k_I{ii,jj} = dN_xi_k_I;
                    obj.dNxi2{ii,jj} = d2N_xi_k_I;
                end
            end
        end
    end
end