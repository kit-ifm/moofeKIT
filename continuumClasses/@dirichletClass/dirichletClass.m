classdef dirichletClass < baseFEClass
    properties (Constant)
        callMassMatrix = false;
        callElements = false;
    end
    properties
        masterObject
        nodeList
        nodalDof = 1;       % 1,2,3,...
        functionDof = [];
        timeFunction = @(t,XYZ) 0;
        time = 0;
        dimension = 3;
    end
    properties (Dependent = true)
        globalNodesDof          % global list of Dirichlet nodes
        qN1
    end

    %% constructor
    methods
        function obj = dirichletClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end
    %% mandatory methods
    methods
%         
        plot(obj,setupObject)
%         
        function initializeShapeFunctions(obj,dofObject)
            % nothing to do here
        end
        function initializeGlobalDofs(obj,dofObject)
            % nothing to do here
        end
        function initializeQV(obj,dofObject)
            % nothing to do here
        end
        function updateGlobalField(obj,dofObject,fieldNameCell)
            for index2 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index2};
                assert(ischar(fieldName));
                if isprop(obj,fieldName)
                    dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
                end
            end
        end
        function updateContinuumFieldPreNewtonLoop(obj,dofObject,fieldName)
            % nothing to do here
        end
        function updateTimeDependentFieldPreNewtonLoop(obj,dofObject,fieldName,time)
            % nothing to do here
        end
        function updateTimeDependentFieldNewtonLoop(obj,dofObject,fieldName,time)
            assert(ischar(fieldName));
            obj.time = time;
            if isprop(obj,fieldName)
                dofObject.(fieldName)(obj.globalNodesDof) = obj.(fieldName);
            end
        end
        function updateContinuumFieldNewtonLoop(obj,dofObject,setupObject,delta,fieldName)
            % nothing to do here
        end
        function updateContinuumFieldPostNewtonLoop(obj,dofObject,setupObject,fieldName)
            % nothing to do here
        end
    end
    %% get und set methods
    methods
        function set.masterObject(obj, input)
            assert(isa(input, 'solidSuperClass'), 'masterObject of Dirichlet Object must be of type solidSuperClass!');
            obj.masterObject = input;
        end
        function set.timeFunction(obj,input)
            assert(isa(input,'function_handle'),'timeFunction must be of type function_handle');
            obj.timeFunction = input;
        end
        function set.nodeList(obj,input)
            checkIntegerMatrix(input);
            obj.nodeList = int32(input(:));
        end
        function set.nodalDof(obj,input)
            checkIntegerMatrix(input);
            obj.nodalDof = input;
        end
        function set.time(obj,t)
            assert(isnumeric(t),'time need to be of class double');
            obj.time = t; 
        end        
        function out = get.globalNodesDof(obj)   
            if ~isempty(obj.masterObject.meshObject.globalNodesDof)
                list = obj.masterObject.meshObject.globalNodesDof(obj.nodeList,obj.nodalDof);
                out = list(:);
            else
                error('globalNodesDof of masterObject is empty. Consider initialization via dofObject.initializeGlobalDofs!');
            end
        end
        function out = get.qN1(obj)
            if numel(obj.nodalDof)==1
                %% FIXME
                if isempty(obj.functionDof)
                    out = obj.timeFunction(obj.time,obj.masterObject.qR(obj.nodeList,obj.nodalDof));
                else
                    out = obj.timeFunction(obj.time,obj.masterObject.qR(obj.nodeList,obj.functionDof));
                end
            else
                out = obj.masterObject.qR(obj.nodeList,obj.nodalDof) + obj.timeFunction(obj.time);
            end
        end
    end
end