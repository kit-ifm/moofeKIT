classdef setupClass < matlab.mixin.Copyable
    properties
        newton = struct('step', [], 'tolerance', 1e-8, 'maximumSteps', 50, 'maximumResiduumNorm', 1e12, 'enforceIteration', false);
        newtonVisco = struct('step', 0, 'tolerance', 1e-8, 'maximumSteps', 50, 'maximumResiduumNorm', 1e12, 'enforceIteration', false);
        toleranceDetJ = 10 * eps;
        usePreconditioning = false;
        totalTime = 1;
        totalTimeSteps = 1;
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
        function set.newton(obj, input)
            % assert(mod(input.step, 1) == 0, "The step specification for Newton's method must be of type integer!");
            assert(isnumeric(input.tolerance) & input.tolerance > 0 & input.tolerance < 1, "The tolerance for Newton's method must be a number between 0 and 1!");
            assert(mod(input.maximumSteps, 1) == 0, "The maximum number of steps for Newton's method must be of type integer!");
            assert(isnumeric(input.maximumResiduumNorm), "The maximum residuum norm for Newton's method must be of type double!");
            obj.newton = input;
        end
        function set.newtonVisco(obj,input)
            assert(mod(input.step, 1) == 0, "The step specification for Newton's method must be of type integer!");
            assert(isnumeric(input.tolerance) & input.tolerance > 0 & input.tolerance < 1, "The tolerance for Newton's method must be a number between 0 and 1!");
            assert(mod(input.maximumSteps, 1) == 0, "The maximum number of steps for Newton's method must be of type integer!");
            assert(isnumeric(input.maximumResiduumNorm), "The maximum residuum norm for Newton's method must be of type double!");
            obj.newtonVisco = input;
        end
        function set.usePreconditioning(obj, val)
            assert(islogical(val), 'usePreconditioning must be logical')
            obj.usePreconditioning = val;
        end
        function set.totalTime(obj, input)
            assert(isnumeric(input) && input > 0, 'totalTime needs to be a positive double!');
            obj.totalTime = input;
        end
        function set.totalTimeSteps(obj, input)
            assert(mod(input, 1) == 0 & input > 0, 'totalTimeSteps needs to be a positive integer!');
            obj.totalTimeSteps = input;
        end
        function set.integrator(obj, input)
            assert(ischar(input), 'Input must be of type string!');
            obj.integrator = input;
        end
        function initialize(obj)
            % initialize newton.step vector size
            obj.newton.step = zeros(obj.totalTimeSteps,1);
            % initialize saveObject
            initialize(obj.saveObject, obj);
        end
        % get methods
        function out = get.timeStepSize(obj)
            % out = round(obj.totalTime/obj.totalTimeSteps, obj.roundDigit);
            out = obj.totalTime/obj.totalTimeSteps;
        end
        function out = get.factorIntegrator(obj)
            if strcmpi(obj.integrator, 'Endpoint')
                out = [1, 1, 0];
            elseif strcmpi(obj.integrator, 'Midpoint') || strcmpi(obj.integrator, 'DiscreteGradient') || strcmpi(obj.integrator, 'LinearImplicit')
                out = [2, 2, 1];
            else
                error('Integrator is not implemented.')
            end
        end
    end
end