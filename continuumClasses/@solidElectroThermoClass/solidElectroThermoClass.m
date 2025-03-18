classdef solidElectroThermoClass < solidSuperClass  
    properties
        additionalFields = 2;
        dofsPerAdditionalField = [1; 1];
    end
    properties (Dependent = true)
        permutationMatrix
    end    
    methods
        %% constructor
        function obj = solidElectroThermoClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
        function out = get.permutationMatrix(obj)
            %% permutation matrix
            % [x1 y1 z1 x2 y2 z2 ... x8 y8 z8 phi1 phi2 phi3 phi4 ... phi8
            %                              24                           32
            % permute to (cf. edof and node vector)
            % [x1 y1 z1 phi1 x2 y2 z2 phi2 ... x8 y8 z8 phi8]
            % index = [1 34 67 101 134 167 201 234 267 301 334 367 401 434 467 501 534 567 601 634 667 701 734 767 772 808 844 880 916 952 988 1024];
            % External dofs
            dimension = obj.dimension;
            numberOfPrimalDOFs = size(obj.meshObject.globalFullEdof,2)-size(obj.mixedFEObject.globalEdof,2);
            mechanicalDOFs = 1:numberOfPrimalDOFs;
            mechanicalDOFs([dimension+1:dimension+2:numberOfPrimalDOFs,dimension+2:dimension+2:numberOfPrimalDOFs]) = [];
            electricalDOFs = dimension+1:dimension+2:numberOfPrimalDOFs;
            temperatureDOFs = dimension+2:dimension+2:numberOfPrimalDOFs;
            index = [mechanicalDOFs, electricalDOFs, temperatureDOFs];
            for i = 2:numel(index)
                index(i) = index(i) + numberOfPrimalDOFs*(i-1);
            end
            out = sparse(numberOfPrimalDOFs,numberOfPrimalDOFs);
            out(index) = 1;
        end
    end
end