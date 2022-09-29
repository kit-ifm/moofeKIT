classdef (Abstract) baseFEClass < handle
    properties
        nodes
        edof                        % local node based edof (element degrees of freedom table)
        globalFullEdof              % global dof based edof
        globalNodesDof              % global node based dof (degrees of freedom)
        dimension
    end
    methods
        %% set and get methods
        function set.edof(obj,input)
            checkIntegerMatrix(input);
            obj.edof = int32(input);
        end
        function set.globalFullEdof(obj,input)
            checkIntegerMatrix(input);
            obj.globalFullEdof = int32(input);
        end
        function set.globalNodesDof(obj,input)
            checkIntegerMatrix(input);
            obj.globalNodesDof = int32(input);
        end
        function set.dimension(obj,input)
            checkIntegerMatrix(input);
            assert((input <= 3) & (input >= 1));
            obj.dimension = input;
        end
    end
end