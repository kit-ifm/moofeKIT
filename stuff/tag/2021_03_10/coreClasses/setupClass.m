classdef setupClass < handle
% handle class
    properties
        fileName = 'save';
        % TODO struct: newton.tolerance; newton.maximumSteps;
        % newton.maximumResiduumNorm
        newtonTolerance = 1e-8;
        newtonMaximumSteps = 50;
        newtonMaximumResiduumNorm = 1e12;
        totalTime = 1;
        totalTimeSteps = 1;
        integrator = 'Endpoint';
        meshType = 'meshGenerator'; %abaqusInput
        meshName = 'Square';
        timeFunction = @(t) obj.timeStepSize
        plotFlag = true;
    end
    properties (Dependent = true)
        timeStepSize
        factorIntegrator
    end
    methods
        % set methods
        function set.fileName(obj,input)
            assert(ischar(input),'Input need to be of type integer!');
            obj.fileName = input;
        end        
        function set.newtonTolerance(obj,input)
            assert(isnumeric(input) & input > 0 & input < 1, 'Input need to a number inbetween 0 and 1!');
            obj.newtonTolerance = input;
        end
        function set.newtonMaximumSteps(obj,input)
            assert(isintger(input),'Input need to be of type integer!');
            obj.newtonMaximumSteps = input;
        end
        function set.newtonMaximumResiduumNorm(obj,input)
            assert(isnumeric(input),'Input need to be of type double!');
            obj.newtonMaximumResiduumNorm = input;
        end
        function set.totalTime(obj,input)
            assert(isnumber(input) & input > 0,'Input need to a postive number!');
            obj.totalTime = input;
        end
        function set.totalTimeSteps(obj,input)
            assert(mod(input,1) == 0,'totalTimeSteps need to be a positiv integer');            
            obj.totalTimeSteps = input;
        end
        function set.integrator(obj,input)
            assert(ischar(input) ,'Input must be of type string!');
            obj.integrator = input;
        end
        % get methods
        function out = get.timeStepSize(obj)
            out = obj.totalTime/obj.totalTimeSteps;
        end
        function out = get.factorIntegrator(obj)
            if strcmpi(obj.integrator,'Endpoint')
                out = [1,0];
            elseif strcmpi(obj.integrator,'Endpoint') || strcmpi(obj.integrator,'DiscreteGradient')
                out = [2,2];
            end
        end
    end
end