classdef solidClass < solidSuperClass  
    properties
        additionalFields = 0;
        dofsPerAdditionalField = [];
    end
    properties
        % additionalFeature
        flagHistoryFields = false;
        historyN % internal variable at timestep N
        historyN1 % internal variable at timestep N1
    end
    methods
        %% constructor
        function obj = solidClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
        function initializeHistoryN1(obj)
            if strcmpi(obj.elementDisplacementType,'displacementSC') && (strcmpi(obj.elementNameAdditionalSpecification,'pHCGJ') || strcmpi(obj.elementNameAdditionalSpecification,'pHCGc') || strcmpi(obj.elementNameAdditionalSpecification,'pHC') || strcmpi(obj.elementNameAdditionalSpecification,'pHCTilde'))
                if isempty(obj.historyN1)
                    numberOfElements = size(obj.meshObject.edof, 1);
                    numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                    for e = 1:numberOfElements
                        for k = 1:numberOfGausspoints
                            obj.historyN1(e,k).C = eye(3);
                            obj.historyN1(e,k).G = eye(3);
                            if (strcmpi(obj.elementNameAdditionalSpecification,'pHCGJ')  || strcmpi(obj.elementNameAdditionalSpecification,'pHCTilde'))
                                obj.historyN1(e,k).J = 1;
                            elseif strcmpi(obj.elementNameAdditionalSpecification,'pHCGc')
                                obj.historyN1(e,k).c = 1;
                            end
                        end
                    end
                end
            end
        end
        function updateHistoryField(obj, dofObject)
            if strcmpi(obj.elementDisplacementType,'displacementSC') && (strcmpi(obj.elementNameAdditionalSpecification,'pHCGJ') || strcmpi(obj.elementNameAdditionalSpecification,'pHCGc') || strcmpi(obj.elementNameAdditionalSpecification,'pHC') || strcmpi(obj.elementNameAdditionalSpecification,'pHCTilde'))
                numberOfElements = size(obj.meshObject.edof, 1);
                numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                for e = 1:numberOfElements
                    for k = 1:numberOfGausspoints
                        obj.historyN(e,k).C = obj.historyN1(e,k).C;
                        obj.historyN(e,k).G = obj.historyN1(e,k).G;
                        if (strcmpi(obj.elementNameAdditionalSpecification,'pHCGJ') || strcmpi(obj.elementNameAdditionalSpecification,'pHCTilde'))
                            obj.historyN(e,k).J = obj.historyN1(e,k).J;
                        elseif strcmpi(obj.elementNameAdditionalSpecification,'pHCGc')
                            obj.historyN(e,k).c = obj.historyN1(e,k).c;
                        end
                    end
                end
            end
        end        
    end
end