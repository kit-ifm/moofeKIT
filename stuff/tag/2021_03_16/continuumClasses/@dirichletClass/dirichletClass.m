classdef dirichletClass < handle
    properties (Constant)
        flagNewton = struct('initializeGlobalDofs',false,...
                            'initializeQV',false,...
                            'initializeShapeFunctions',false,...
                            'updateTimeDependentField',true,...
                            'updateContinuumFieldPreNewtonLoop',false,...
                            'callElements',false,...
                            'updateContinuumFieldNewtonLoop',false);
    end
    properties
        masterObject
        nodeList
        nodalDof = 1;       % 1,2,3         
        timeFunction = @(t,XYZ) 0;
        time = 0;
    end
    properties (Dependent = true)
        globalDirichletList
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
    
    %% get und set methods
    methods
        function set.masterObject(obj,input)
            assert(isa(input,'solidClass')||isa(input,'solidThermoClass'),'masterObject must be compatibel to dirichlet class');
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
        function out = get.globalDirichletList(obj)   
            if ~isempty(obj.masterObject.globalNodesDof)
                list = obj.masterObject.globalNodesDof(obj.nodeList,obj.nodalDof);
                out = list(:);
            else
                error('globalNodesDof of masterObject is empty. Consider initialization via dofObject.initializeGlobalDofs!');
            end
        end
        function out = get.qN1(obj)
            if numel(obj.nodalDof)==1
                out = obj.timeFunction(obj.time,obj.masterObject.qR(obj.nodeList,obj.nodalDof));
            else
                out = obj.masterObject.qR(obj.nodeList,obj.nodalDof) + obj.timeFunction(obj.time);
            end
        end
    end
end