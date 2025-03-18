classdef axisymmetricSolidClass < solidSuperClass  
    properties
        additionalFields = 0;
        dofsPerAdditionalField = [];
    end
    methods
        %% constructor
        function obj = axisymmetricSolidClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end
end