classdef stringClass < solidSuperClass
    properties (Constant = true)
        additionalFields = 0;
    end
    properties (SetAccess = private)
        numberOfDofsPerNode
    end
    methods

        %% constructor
        function obj = stringClass(dofObject,numberOfSpatialDimensions)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            updateNumberOfDofsPerNode(obj,numberOfSpatialDimensions)
        end
        
        function initializeGlobalDofs(obj, dofObject)
            numberOfDofPerNode = size(obj.meshObject.nodes,2);
            edof = obj.meshObject.edof;
            numberOfNodes = size(obj.meshObject.nodes,1);
            numberOfElements = size(edof,1);
            globalFullEdof = zeros(numberOfElements,size(edof,2)*numberOfDofPerNode);
            for j = 1:numberOfDofPerNode
                globalFullEdof(:,j:numberOfDofPerNode:end) = dofObject.totalNumberOfDofs + edof*numberOfDofPerNode-(numberOfDofPerNode-j);
            end
            globalNodesDof = zeros(numberOfNodes,numberOfDofPerNode);
            globalNodesDof(:,1:numberOfDofPerNode) = dofObject.totalNumberOfDofs + kron(ones(numberOfNodes,1),1:numberOfDofPerNode)+kron((0:numberOfDofPerNode:numberOfDofPerNode*numberOfNodes-1)',ones(1,numberOfDofPerNode));
            obj.meshObject.globalFullEdof = globalFullEdof;
            obj.meshObject.globalNodesDof = globalNodesDof;
            dofObject.totalNumberOfDofs = dofObject.totalNumberOfDofs + numberOfNodes*numberOfDofPerNode;
            % mixed elements
            if ~contains(obj.elementDisplacementType,{'displacement'})
                obj.mixedFEObject.initializeMixedElements(dofObject,obj);
                obj.meshObject.globalFullEdof = [obj.meshObject.globalFullEdof, obj.mixedFEObject.globalEdof];
            end
        end

        
    end
    methods
        function set.numberOfDofsPerNode(obj, input)
            assert(mod(input, 1) == 0, 'Input must be of type integer!')
            obj.numberOfDofsPerNode = input;
        end
        function value = get.numberOfDofsPerNode(obj)
            value = obj.numberOfDofsPerNode;
        end
    end
    methods (Access=private)
        function updateNumberOfDofsPerNode(obj,nDIM)           
            obj.numberOfDofsPerNode = nDIM;
        end
    end
end