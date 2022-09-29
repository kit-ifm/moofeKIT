classdef solidSuperClass < baseFEClass
    properties (Abstract = true, Constant = true)
        flagNewton
        additionalFields
    end
    properties
        % mesh
        orderShapeFunctions
        numberOfGausspoints
        shapeFunctions
        % material
        materialName
        materialData
        % solver
        qN
        qN1
        vN
        vN1
    end
    properties (Dependent = true)
        qR
        massMatrix
    end    
    methods
        function out = get.qR(obj)
            out = obj.nodes;
        end
        function out = get.massMatrix(obj)
            out = massMatrixElement(obj);
        end
    end
end