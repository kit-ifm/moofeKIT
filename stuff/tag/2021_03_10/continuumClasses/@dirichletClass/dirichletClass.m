classdef dirichletClass < handle
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
    properties (Constant)
        flagNewton = struct('initializeGlobalDofs',false,...
                            'initializeQV',false,...
                            'initializeShapeFunctions',false,...
                            'updateTimeDependentField',true,...
                            'updateContinuumFieldPreNewtonLoop',false,...
                            'callElements',false,...
                            'updateContinuumFieldNewtonLoop',false);
    end
    
    %% constructor
    methods
        function obj = dirichletClass
        end
    end
    
    %% get und set methods
    methods
        function set.masterObject(obj,input)
            assert(isa(input,'solidClass'),'masterObject must be compatibel to dirichlet class');
            obj.masterObject = input;
        end
        function set.timeFunction(obj,input)
%             assert(isa(input,'function_handle'),'timeFunction must be of type function_handle');
            obj.timeFunction = input;
        end
        function set.nodeList(obj,input)
            checkIntegerMatrix(input);
            obj.nodeList = int32(input(:));
        end
        function set.nodalDof(obj,input)
            checkIntegerMatrix(input);
%             assert((input <= 3) & (input >= 1));
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
                error('globalNodesDof of masterObject is empty. Consider initialization via dofsObject.initializeGlobalDofs!');
            end
        end
        function out = get.qN1(obj)
            if numel(obj.nodalDof)==1
                out = obj.timeFunction(obj.time,obj.masterObject.qR(obj.nodeList,obj.nodalDof));
%                 TODO add other cases (2:3, 1:3) etc.
            else
                out = obj.masterObject.qR(obj.nodeList,obj.nodalDof) + obj.timeFunction(obj.time);
            end
        end
%         function updateField(obj,dofsObject,fieldNameCell,flagUpdate)
%             for index1 = 1:numel(fieldNameCell)
%                 fieldName = fieldNameCell{index1};
%                 field = dofsObject.(fieldName);
%                 assert(ischar(fieldName) & ischar(flagUpdate));
%                 for index2 = 1:numel(obj)
%                     if  isfield(obj(index2).flagNewton,flagUpdate)
%                         if obj(index2).flagNewton.(flagUpdate)
%                             obj(index2).(fieldName) = field(obj(index2).globalNodesDof);
%                         end
%                     end
%                 end
%             end            
%             assert(ischar(fieldName) & isnumeric(field) & ischar(flagUpdate));
%             if isprop(obj,fieldName)
%                 for index = 1:numel(obj)
%                     if  isfield(obj(index).flagNewton,flagUpdate)
%                         if obj(index).flagNewton.(flagUpdate)
%                             obj(index).(fieldName) = field(obj(index).globalNodesDof);
%                         end
%                     end
%                 end
%             end
%         end
    end
end