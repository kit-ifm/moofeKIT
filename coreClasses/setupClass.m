classdef setupClass < matlab.mixin.Copyable
    properties
        newton = struct('step', 0, 'tolerance', 1e-8, 'maximumSteps', 50, 'maximumResiduumNorm', 1e12, 'enforceIteration', false);
        toleranceDetJ = 10*eps;
        usePreconditioning = false;
        totalTime = 1;
        totalTimeSteps = 1;
        roundDigit = 4;
        integrator = 'Endpoint';
        meshType = 'meshGenerator'; %abaqusInput
        meshName = 'Square';
        timeFunction = @(t) obj.timeStepSize
        timeStep
        time = 0;
        % load objects
        saveObject
        plotObject
        % 
        computePostData = false;
    end
    properties (Dependent = true)
        timeStepSize
        factorIntegrator
    end
    methods
        %% constructor
        function obj = setupClass()
            obj.saveObject = saveClass;
            obj.plotObject = plotClass;
        end
        % set methods
        function set.newton(obj,input)
            assert(mod(input.step, 1) == 0, "The step specification for Newton's method must be of type integer!");
            assert(isnumeric(input.tolerance) & input.tolerance > 0 & input.tolerance < 1, "The tolerance for Newton's method must be a number between 0 and 1!");
            assert(mod(input.maximumSteps, 1) == 0, "The maximum number of steps for Newton's method must be of type integer!");
            assert(isnumeric(input.maximumResiduumNorm), "The maximum residuum norm for Newton's method must be of type double!");
            obj.newton = input;
        end
        function set.usePreconditioning(obj, val)
            assert(islogical(val), 'usePreconditioning must be logical')
            obj.usePreconditioning = val;
        end
        function set.totalTime(obj,input)
            assert(isnumeric(input) & input > 0,'Input need to a postive number!');
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
            out = round(obj.totalTime/obj.totalTimeSteps,obj.roundDigit);
        end
        function out = get.factorIntegrator(obj)
            if strcmpi(obj.integrator,'Endpoint')
                out = [1,1,0];
            elseif strcmpi(obj.integrator,'Midpoint') || strcmpi(obj.integrator,'DiscreteGradient')
                out = [2,2,1];
            else
                error('Integrator is not implemented.')
            end
        end
    end
end