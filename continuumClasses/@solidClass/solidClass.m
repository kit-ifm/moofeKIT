classdef solidClass < solidSuperClass  
    properties (Constant = true)
        additionalFields = 0;
    end
    methods
        %% constructor
        function obj = solidClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end
end