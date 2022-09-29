classdef numericalTangentClass < matlab.mixin.Copyable
    properties
        computeNumericalTangent = false;
        showDifferences = false;
        type = 'standard';
        epsilon = 1e-6;
        flag = false;
        flagDiscreteGradient
    end
    methods

        %% set and get methods
        function set.computeNumericalTangent(obj, val)
            assert(islogical(val), 'computeNumericalTangent must be logical')
            obj.computeNumericalTangent = val;
        end

        function set.showDifferences(obj, val)
            assert(islogical(val), 'showDifferences must be logical')
            obj.showDifferences = val;
        end

        function set.type(obj, val)
            assert(ischar(val), 'type must be of type string!');
            obj.type = val;

            updateEpsilon(obj, val);
        end

        function set.epsilon(obj, val)
            assert(isnumeric(val), 'epsilon must be of type double!');
            obj.epsilon = val;
        end

        function set.flag(obj, val)
            assert(islogical(val), 'flag must be logical')
            obj.flag = val;
        end

        function updateEpsilon(obj, type)
            if strcmpi(type, 'standard')
                obj.epsilon = 1e-6;
            elseif strcmpi(type, 'centralDifferences')
                obj.epsilon = sqrt(eps);
            elseif strcmpi(type, 'complex')
                obj.epsilon = 1i * 1e-50;
            end
        end

        function initializeFlagDiscreteGradient(obj, numberOfGausspoints, numberOfDiscreteDerivatives)
            assert(mod(numberOfGausspoints, 1) == 0, 'numberOfGausspoints must be an integer');
            assert(mod(numberOfDiscreteDerivatives, 1) == 0, 'numberOfDiscreteDerivatives must be an integer');
            obj.flagDiscreteGradient = false(numberOfGausspoints, numberOfDiscreteDerivatives);
        end
    end
end
