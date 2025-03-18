classdef beamClass < solidSuperClass
    properties
        additionalFields = 1;
        dofsPerAdditionalField = 1;
    end
    properties
        theory  % beam theory: 'bernoulli', 'timoshenko' or 'geometricallyExact'
        selectiveReducedShapeFunctionObject
        historyVariableN % internal strain variable at timestep N
        historyVariableN1 % internal strain variable at timestep N+1
    end
    properties (Dependent = true)
        permutationMatrix
    end
    properties (SetAccess = private)
        numberOfDofsPerNode = 1;
    end
    methods
        
        %% constructor
        function obj = beamClass(dofObject,varargin)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            if ~isempty(varargin)
                updateNumberOfDisplacementDofsPerNode(obj,varargin{1})
            end
            obj.selectiveReducedShapeFunctionObject = shapeFunctionClass();
        end
        
        function initializeQR(obj)
            % method initializeQR
            numberOfNodes = size(obj.meshObject.nodes, 1);
            if strcmpi(obj.theory , 'GeometricallyExact')
                obj.qR = zeros(numberOfNodes, 3);
                obj.qR(:,1:2) = obj.meshObject.nodes; % 2D position
                obj.qR(:,3) = atan2((obj.qR(end,2)-obj.qR(1,2)),(obj.qR(end,1)-obj.qR(1,1))); % rotation angle (assuming straight config)
            elseif strcmpi(obj.theory , 'Bernoulli') || strcmpi(obj.theory , 'Timoshenko')
                obj.qR = zeros(numberOfNodes, 2); % vertical displacements and rotation angle
            end
        end
        
        function out = getAngularMomentum(obj)   % Get method for the angular momentum in 2D.
            out = 0;
            if strcmpi(obj.theory , 'GeometricallyExact')
                numberOfNodes = size(obj.meshObject.nodes, 1);
                skew_matrix = [0 , 1;
                    -1, 0];
                for j = 1:numberOfNodes
                    out = out + obj.qN1(j,1:2) * skew_matrix * obj.L(j,1:2)' + obj.L(j,3); % first term: r \cross rhoA*v, second term: rhoI*omega
                end
            end
        end
        
        function out = get.permutationMatrix(obj)
            
            %% permutation matrix
            % [x1 x2 phi1 phi2]
            % permute to (cf. edof and node vector)
            % [x1 phi1 x2 phi2]
            %
            % or [x1 y2 x2 y2 phi1 phi2]
            % to [x1 y1 phi1 x2 y2 phi2]
            
            % total nbr of DOFs (including mixed fields)
            numberOfPrimalDOFs = size(obj.meshObject.globalFullEdof, 2);
            if ~isempty(obj.mixedFEObject.numberOfDofs)
                % remove dimension of mixed quantities from this matrix
                numberOfPrimalDOFs = numberOfPrimalDOFs - obj.mixedFEObject.numberOfDofs;
            end
            % get DOFs related to rotation
            nbr_of_translation_dofs = obj.numberOfDofsPerNode;
            phiDOFs = nbr_of_translation_dofs+1:nbr_of_translation_dofs + 1:numberOfPrimalDOFs;
            % get DOFs related to translation
            displacementDOFs = 1:numberOfPrimalDOFs;
            displacementDOFs(nbr_of_translation_dofs+1:nbr_of_translation_dofs+1:numberOfPrimalDOFs) = [];
            index = [displacementDOFs, phiDOFs];
            for i = 2:numel(index)
                index(i) = index(i) + numberOfPrimalDOFs * (i - 1);
            end
            out = sparse(numberOfPrimalDOFs, numberOfPrimalDOFs);
            out(index) = 1;
        end
        
        function set.theory(obj, input)
            assert(ischar(input), 'Input must be of type string!');
            if strcmp(input, 'bernoulli')
                input = 'Bernoulli';
            elseif strcmp(input, 'timoshenko')
                input = 'Timoshenko';
            elseif strcmp(input, 'geometricallyExact')
                input = 'GeometricallyExact';
            end
            if any(strcmp(input, {'Bernoulli', 'Timoshenko', 'GeometricallyExact'}))
                obj.theory = input;
            else
                error('This beam theory is not implemented yet!');
            end
        end

        function initializeHistoryVariableN1(obj)
            if isempty(obj.historyVariableN1) && strcmp(obj.elementNameAdditionalSpecification,'PhIrreducible')
                numberOfElements = size(obj.meshObject.edof, 1);
                numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                obj.historyVariableN1 = zeros(numberOfElements, numberOfGausspoints, 2); % strain w.r.t elongation and shear is initially zero
            end
            if isempty(obj.historyVariableN) && strcmp(obj.elementNameAdditionalSpecification,'PhIrreducible')
                obj.historyVariableN = obj.historyVariableN1; % strain w.r.t elongation and shear is initially zero
            end
        end
        
        function updateIteratedHistoryField(obj,setupObject)
            [~] = callElements(obj, setupObject, 'postData');
        end

        function updateHistoryField(obj, ~)
            obj.historyVariableN = obj.historyVariableN1;
        end
        
    end
    methods (Access=private)
        function updateNumberOfDisplacementDofsPerNode(obj,nDOF)
            obj.numberOfDofsPerNode = nDOF;
        end
    end
end