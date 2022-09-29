classdef dirichletClass < handle
    properties
        masterObject
        masterNodesDofDirichlet
        timeFunction = @(t) 0;
        time = 0;
        masterNodesDirection = 1;       % 1,2,3 
    end
    properties (Dependent = true)
        masterGlobalNodesDof
        qN1
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
        function set.masterNodesDofDirichlet(obj,input)
            checkIntegerMatrix(input);
            obj.masterNodesDofDirichlet = int32(input(:));
        end
        function set.masterNodesDirection(obj,input)
            checkIntegerMatrix(input);
            assert((input <= 3) & (input >= 1));
            obj.masterNodesDirection = input;
        end
        function set.time(obj,t)
            assert(isnumeric(t),'time need to be of class double');
            obj.time = t; 
        end        
        function out = get.masterGlobalNodesDof(obj)
            countBC = 0;
            for ii=1:numel(obj)
                countBC = countBC + numel(obj(ii).masterNodesDofDirichlet);
            end                
            out = zeros(countBC,1);
            countBC = 0;
            for ii=1:numel(obj)
                dofs = obj(ii).masterObject.globalNodesDOF(obj(ii).masterNodesDofDirichlet,obj(ii).masterNodesDirection);
                out(countBC+1:countBC+size(dofs,1)) = dofs;
                countBC = countBC + size(dofs,1);
            end
        end
        function out = get.qN1(obj)
            out = zeros(size(obj.masterNodesDofDirichlet,1),1);
            out(:) = obj.masterObject.QR(obj.masterNodesDofDirichlet,obj.masterNodesDirection) + obj.timeFunction(obj.time);
        end     
    end
end