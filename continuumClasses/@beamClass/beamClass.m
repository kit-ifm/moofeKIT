classdef beamClass < solidSuperClass
    properties (Constant = true)
        additionalFields = 1;
    end
    properties (Dependent = true)
        permutationMatrix
    end
    methods

        %% constructor
        function obj = beamClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end

        function initializeQR(obj)
            % method initializeQR
            numberOfNodes = size(obj.meshObject.nodes, 1);
            obj.qR = zeros(numberOfNodes, 2);
        end

        function out = get.permutationMatrix(obj)

            %% permutation matrix
            % [x1 x2 phi1 phi2]
            % permute to (cf. edof and node vector)
            % [x1 phi1 x2 phi2]
            % External dofs
            dimension = obj.dimension;
            numberOfPrimalDOFs = size(obj.meshObject.globalFullEdof, 2);
            phiDOFs = dimension+1:dimension + 1:numberOfPrimalDOFs;
            displacementDOFs = 1:numberOfPrimalDOFs;
            displacementDOFs(dimension+1:dimension+1:numberOfPrimalDOFs) = [];
            index = [displacementDOFs, phiDOFs];
            for i = 2:numel(index)
                index(i) = index(i) + numberOfPrimalDOFs * (i - 1);
            end
            out = sparse(numberOfPrimalDOFs, numberOfPrimalDOFs);
            out(index) = 1;
        end
    end
end