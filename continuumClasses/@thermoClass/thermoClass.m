classdef thermoClass < solidSuperClass
    properties
        additionalFields = 0;
        dofsPerAdditionalField = [];
    end
    methods

        %% constructor
        function obj = thermoClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end

        function initializeQR(obj)
            % method initializeQR
            numberOfNodes = size(obj.meshObject.nodes, 1);
            obj.qR = zeros(numberOfNodes, 1);
        end
    end
end