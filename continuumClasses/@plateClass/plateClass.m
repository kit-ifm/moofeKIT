classdef plateClass < solidSuperClass
    properties
        additionalFields = 0;
        dofsPerAdditionalField = [];
    end
    properties
        theory = 'mindlin'  % plate theory: 'mindlin' or 'kirchhoff'
        h = 1;              % plate thickness
    end
    properties (SetAccess = private)
        numberOfDofsPerNode
    end
    methods

        %% constructor
        function obj = plateClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            
            updateNumberOfDofsPerNode(obj)
        end
        
        function initializeGlobalDofs(obj, dofObject)
            numberOfDofPerNode = obj.numberOfDofsPerNode;
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
            if ~contains(obj.elementDisplacementType,{'displacement', 'thermo', 'beam'})
                obj.mixedFEObject.initializeMixedElements(dofObject,obj);
                obj.meshObject.globalFullEdof = [obj.meshObject.globalFullEdof, obj.mixedFEObject.globalEdof];
            end
        end

        function initializeQR(obj)
            % method initializeQR
            numberOfNodes = size(obj.meshObject.nodes, 1);
            obj.qR = zeros(numberOfNodes, obj.numberOfDofsPerNode);
        end
        
    end
    methods
        function set.h(obj, input)
            assert(isnumeric(input), 'Plate thickness must be of type float!');
            assert(input > 0, 'Plate thickness must be larger than zero!');
            obj.h = input;
        end
        function value = get.h(obj)
            value = obj.h;
        end
        function set.theory(obj, input)
            assert(ischar(input), 'Input must be of type string!');
            input = lower(input);
            if any(strcmp(input, {'mindlin', 'kirchhoff'}))
                obj.theory = input;
                updateNumberOfDofsPerNode(obj);
            else
                error('This plate theory is not implemented yet!');
            end
        end
        function set.numberOfDofsPerNode(obj, input)
            assert(mod(input, 1) == 0, 'Input must be of type integer!')
            obj.numberOfDofsPerNode = input;
        end
        function value = get.numberOfDofsPerNode(obj)
            value = obj.numberOfDofsPerNode;
        end
    end
    methods (Access=private)
        function updateNumberOfDofsPerNode(obj)
            if strcmpi(obj.theory, 'mindlin')
                obj.numberOfDofsPerNode = 3;
            elseif strcmpi(obj.theory, 'kirchhoff')
                obj.numberOfDofsPerNode = 1;
            end
        end
    end
end