classdef solidVelocityClass < solidSuperClass
    properties (Dependent = true)
        additionalFields
        dofsPerAdditionalField
        permutationMatrix
    end
    properties
        % additionalFeature
        flagHistoryFields = false;
        historyN % internal variable at timestep N
        historyN1 % internal variable at timestep N1
    end
    methods

        %% constructor
        function obj = solidVelocityClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end

        function out = get.additionalFields(obj)
            out = 1;
        end

        function out = get.dofsPerAdditionalField(obj)
            out = obj.dimension;
        end

        function out = get.permutationMatrix(obj)
            %% permutation matrix
            % [x1 y1 z1 x2 y2 z2 ... x8 y8 z8 theta1 theta2 theta3 theta4 ... theta8
            %                              24                           32
            % permute to (cf. edof and node vector)
            % [x1 y1 z1 theta1 x2 y2 z2 theta2 ... x8 y8 z8 theta8]
            % index = [1 34 67 101 134 167 201 234 267 301 334 367 401 434 467 501 534 567 601 634 667 701 734 767 772 808 844 880 916 952 988 1024];
            % External dofs
            dimension = obj.dimension;
            numberOfPrimalDOFs = size(obj.meshObject.globalFullEdof, 2) - size(obj.mixedFEObject.globalEdof, 2);
            displacementDOFs = 1:numberOfPrimalDOFs;
            velocityDOFs = 1:numberOfPrimalDOFs;
            indicesDisplacementDOFS = [];
            indicesVelocityDOFS = [];
            for d=1:dimension
                indicesDisplacementDOFS = [indicesDisplacementDOFS, d:dimension*2:numberOfPrimalDOFs];
                indicesVelocityDOFS = [indicesVelocityDOFS, dimension+d:dimension*2:numberOfPrimalDOFs];
            end
            displacementDOFs(indicesVelocityDOFS) = [];
            velocityDOFs(indicesDisplacementDOFS) = [];

            index = [displacementDOFs, velocityDOFs];
            for i = 2:numel(index)
                index(i) = index(i) + numberOfPrimalDOFs * (i - 1);
            end
            out = sparse(numberOfPrimalDOFs, numberOfPrimalDOFs);
            out(index) = 1;
        end

        function updateContinuumFieldPostNewtonLoop(obj,dofObject,setupObject,fieldName)
        % method                    updateContinuumFieldPostNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
            field = dofObject.qN1;
            assert(ischar(fieldName));
            if isprop(obj,fieldName)
                obj.(fieldName) = field(obj.meshObject.globalNodesDof);
            end
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