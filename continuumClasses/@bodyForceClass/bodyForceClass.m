classdef bodyForceClass < baseFEClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        callMassMatrix = false;
        callElements = true;
    end
    properties (SetAccess = private)
        elementGeometryType
    end
    properties
        time = 0;
        dimension = 3;
        %
        storageFEObject
        shapeFunctionObject% = lagrangeShapeFunctionClass();
        %
        masterObject
        typeOfLoad
        timeFunction
        loadFunction
        elementData = struct();
    end
    methods
        %% Constructor
        function obj = bodyForceClass(dofObject)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            obj.storageFEObject = storageFEClass();
            obj.shapeFunctionObject = shapeFunctionClass();
            obj.meshObject = meshClass();
        end

        %% FE- Routines
        out = deadLoad(obj,varargin)

        %% Tools
        obj = constructor(obj,mesh,varargin)                 % Constructor.
        plot(obj,varargin)
    end
    %% mandatory methods
    methods
        function initializeShapeFunctions(obj,dofObject)
            shapeFunctionObject = obj.shapeFunctionObject;
            dimension = obj.masterObject.dimension;
            numberOfNodes = size(obj.meshObject.edof,2);
            obj.elementGeometryType = determineElementGeometryType(dimension, numberOfNodes);
            shapeFunctionObject.computeShapeFunction(dimension, numberOfNodes, obj.elementGeometryType);
        end
        function initializeGlobalDofs(obj,dofObject)
            masterObject = obj.masterObject;
            dimension = obj.dimension;
            edof = obj.meshObject.edof;
            numberOfElements = size(edof,1);
            globalFullEdof = zeros(numberOfElements,size(edof,2)*dimension);
            for index2 = 1:dimension
                temp = zeros(numberOfElements,size(edof,2));
                temp(:) = masterObject.meshObject.globalNodesDof(edof,index2);
                globalFullEdof(:,index2:dimension:end) = temp;
            end
            obj.meshObject.globalFullEdof = globalFullEdof;
        end
        function initializeQV(obj,dofObject)
            % nothing to do here
        end
        function updateGlobalField(obj,dofObject,fieldNameCell)
            % nothing to do here
            %
            %             for index2 = 1:numel(fieldNameCell)
            %                 fieldName = fieldNameCell{index2};
            %                 assert(ischar(fieldName));
            %                 if isprop(obj,fieldName)
            %                     dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
            %                 end
            %             end
        end
        function updateContinuumFieldPreNewtonLoop(obj,setupObject,fieldNameCell)
            % nothing to do here
        end
        function updateTimeDependentFieldPreNewtonLoop(obj,dofObject,fieldName,time)
            assert(ischar(fieldName));
            obj.time = time;
            if isprop(obj,fieldName)
                dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
            end
        end
        function updateTimeDependentFieldNewtonLoop(obj,dofObject,fieldName,time)
            % nothing to do here
        end
        function updateContinuumFieldNewtonLoop(obj,dofObject,setupObject,delta,fieldName)
            % nothing to do here
        end
        function updateContinuumFieldPostNewtonLoop(obj,dofObject,setupObject,fieldName)
            % nothing to do here
        end
    end
end
