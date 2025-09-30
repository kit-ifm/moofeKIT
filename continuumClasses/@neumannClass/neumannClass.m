classdef neumannClass < baseFEClass
    properties (Constant)
        callMassMatrix = false;
        callElements = true;
    end
    properties
        storageFEObject
        shapeFunctionObject
        masterObject % object of solidSuperClass
        loadVector % load vector (defined in the global coordinate system)
        loadVectorFunction = @(XYZ) 0; % load vector as a function of the global coordinates (can be used alternatively to constant 'loadVector')
        loadGeometry % 'line' or 'area'
        loadType = 'deadLoad'; % 'deadLoad' or 'followerLoad'
        loadPhysics = 'mechanical'; % 'mechanical', 'thermal', 'electrical', tbc
        timeFunction = @(t) 1;
        projectionType = 'none' % which length is used for calculation of area (not needed in 1D) --> 'none': real length, 'x': length in xdirection, 'y': length in ydirection

        elementData = struct('externalEnergy', 0);
    end
    properties (SetAccess = private)
        elementGeometryType
        masterNodes
        nodalForces
        globalEdof
    end
    properties (Dependent = true)
        numberOfElements
    end

    %% constructor
    methods
        function obj = neumannClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            obj.storageFEObject = storageFEClass();
            obj.shapeFunctionObject = shapeFunctionClass();
        end
    end

    %% mandatory methods
    methods
        function initializeShapeFunctions(obj, dofObject)
            % initialize shape functions
            if strcmp(obj.loadGeometry, 'line')
                dimensionShapeFunctions = 1;
            elseif strcmp(obj.loadGeometry, 'area')
                dimensionShapeFunctions = 2;
            end
            numberOfNodes = size(obj.meshObject.edof, 2);
            if isempty(obj.shapeFunctionObject.numberOfGausspoints)
                obj.shapeFunctionObject.numberOfGausspoints = nthroot(obj.masterObject.shapeFunctionObject.numberOfGausspoints, obj.masterObject.dimension)^dimensionShapeFunctions;
            end
            if isempty(obj.shapeFunctionObject.order)
                obj.shapeFunctionObject.order = obj.masterObject.shapeFunctionObject.order;
            end
            obj.elementGeometryType = determineElementGeometryType(dimensionShapeFunctions, numberOfNodes);
            obj.shapeFunctionObject.computeShapeFunction(dimensionShapeFunctions, numberOfNodes, obj.elementGeometryType);
        end
        function initializeGlobalDofs(obj, dofObject)
            edof = obj.meshObject.edof;
            numberOfAdditionalDofsPerNode = sum(obj.masterObject.dofsPerAdditionalField);
            displacementDofsPerNode = size(obj.masterObject.qR, 2) - numberOfAdditionalDofsPerNode;
            numberOfNodesPerElement = size(edof, 2);
            if strcmp(obj.loadPhysics, 'mechanical')
                numberOfDofsPerElement = numberOfNodesPerElement * displacementDofsPerNode;
                globalEdof = zeros(obj.numberOfElements, numberOfDofsPerElement);
                for ii = 1:displacementDofsPerNode
                    temp = zeros(obj.numberOfElements, numberOfNodesPerElement);
                    temp(:) = obj.masterObject.meshObject.globalNodesDof(edof, ii);
                    globalEdof(:, ii:displacementDofsPerNode:end) = temp;
                end
            elseif strcmp(obj.loadPhysics, 'thermal')
                globalEdof = zeros(obj.numberOfElements, numberOfNodesPerElement);
                if isa(obj.masterObject, 'solidThermoClass')
                    globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, displacementDofsPerNode+1);
                elseif isa(obj.masterObject, 'solidElectroThermoClass')
                    globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, displacementDofsPerNode+2);
                else
                    error('Not implemented yet!');
                end
            elseif strcmp(obj.loadPhysics, 'electrical')
                globalEdof = zeros(obj.numberOfElements, numberOfNodesPerElement);
                if isa(obj.masterObject, 'solidElectroClass') || isa(obj.masterObject, 'solidElectroThermoClass')
                    globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, displacementDofsPerNode+1);
                else
                    error('Not implemented yet!');
                end
            end
            obj.globalEdof = globalEdof;

            determineMasterNodes(obj);

            % verification of input data
            verificationOfInputData(obj)
        end
        function initializeQV(obj, dofObject)
            % nothing to do here
        end
        function updateGlobalField(obj, dofObject, fieldNameCell)
            for index2 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index2};
                assert(ischar(fieldName));
                if isprop(obj, fieldName)
                    dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
                end
            end
        end
        function updateContinuumFieldPreNewtonLoop(obj, setupObject, fieldNameCell)
            % nothing to do here
        end
        function updateTimeDependentFieldPreNewtonLoop(obj, dofObject, fieldName, time)
            assert(ischar(fieldName));
            if isprop(obj, fieldName)
                dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
            end
        end
        function updateTimeDependentFieldNewtonLoop(obj, dofObject, fieldName, time)
            % nothing to do here
        end
        function updateContinuumFieldNewtonLoop(obj, dofObject, delta, setupObject, fieldName)
            % nothing to do here
        end
        function updateContinuumFieldPostNewtonLoop(obj, dofObject, setupObject, fieldName)
            % nothing to do here
        end
    end

    %% get und set methods
    methods

        %% set routines
        function set.masterObject(obj, value)
            assert(isa(value, 'solidSuperClass'), 'masterObject of Neumann Object must be of type solidSuperClass!')
            obj.masterObject = value;
        end

        function set.loadVector(obj, value)
            assert(size(value, 1) == 1 || size(value, 2) == 1, 'loadVector must be a vector!');
            if size(value, 1) == 1
                value = value.';
            end
            obj.loadVector = value;
        end

        function set.loadVectorFunction(obj, value)
            assert(isa(value,'function_handle'), 'loadVectorFunction must be of type function_handle');
            obj.loadVectorFunction = value;
        end

        function set.loadGeometry(obj, value)
            assert(ischar(value), 'loadGeometry must be of type string!')
            value = lower(value);
            if any(strcmp(value, {'line', 'area'}))
                obj.loadGeometry = value;
            else
                error('Not implemented yet for given loadGeometry!');
            end
        end
        function out = get.loadGeometry(obj)
            if isempty(obj.loadGeometry)
                error('loadGeometry is not specified yet!');
            else
                out = obj.loadGeometry;
            end
        end

        function set.loadType(obj, value)
            assert(ischar(value), 'loadType must be of type string!')
            assert(any(strcmp(value, {'deadLoad', 'followerLoad'})), 'loadType must be either "deadLoad" or "followerLoad"!');
            assert(~strcmp(value, 'followerLoad'), 'followerLoad loads are not implemented yet!');
            obj.loadType = value;
        end

        function set.loadPhysics(obj, value)
            assert(ischar(value), 'loadPhysics must be of type string!')
            assert(any(strcmp(value, {'mechanical', 'thermal', 'physical'})), 'loadPhysics must be either "mechanical", "thermal" or "electrical"!');
            obj.loadPhysics = value;
        end

        function set.timeFunction(obj, input)
            assert(isa(input, 'function_handle'), 'timeFunction must be of type function_handle');
            obj.timeFunction = input;
        end

        function set.projectionType(obj, projection)
            assert(any(strcmp(projection, {'none', 'x', 'y', 'z'})), "projectionType must be 'none', 'x', 'y' or 'z'");
            obj.projectionType = projection;
        end


        function out = get.numberOfElements(obj)
            out = size(obj.meshObject.edof, 1);
        end
    end

    %% Further methods
    methods
        function determineMasterNodes(obj)
            obj.masterNodes = unique(obj.meshObject.edof');
        end

        function updateNodalForces(obj, dofObject)
            R = full(sparse(vertcat(obj.storageFEObject.dataFE(:).indexReI), 1, vertcat(obj.storageFEObject.dataFE(:).Re), dofObject.totalNumberOfDofs, 1));
            obj.nodalForces(:, :) = -R(obj.masterObject.meshObject.globalNodesDof(obj.masterNodes, :));
        end

        function verificationOfInputData(obj)
            % load geometry
            if obj.masterObject.dimension == 2 && isa(obj.masterObject, 'solidSuperClass')
                if ~isa(obj.masterObject, 'plateClass')
                    assert(~strcmp(obj.loadGeometry, 'area'), 'Geometry of Neumann object can not be "area" in this case!');
                end
            end

            % load vector
            if strcmp(obj.loadPhysics, 'mechanical')
                % if isa(obj.masterObject,'solidVelocityClass')
                %     displacementDofsPerNode = size(obj.masterObject.qR, 2);
                % else
                    numberOfAdditionalDofsPerNode = sum(obj.masterObject.dofsPerAdditionalField);
                    displacementDofsPerNode = size(obj.masterObject.qR, 2) - numberOfAdditionalDofsPerNode;
                % end
                if ~isempty(obj.loadVector)
                    % check loadVector
                    assert(size(obj.loadVector, 1) == displacementDofsPerNode, ['loadVector must be of size "', num2str(displacementDofsPerNode), ' x 1"']);
                else
                    % check loadVectorFunction
                    assert(size(obj.loadVectorFunction([1;1;1]), 1) == displacementDofsPerNode, ['loadVectorFunction must be of size "', num2str(displacementDofsPerNode), ' x 1"']);
                end
            elseif any(strcmp(obj.loadPhysics, {'thermal', 'electrical'}))
                if ~isempty(obj.loadVector)
                    % check loadVector
                    assert(size(obj.loadVector, 1) == 1, 'loadVector must be of size "1 x 1"');
                else
                    % check loadVectorFunction
                    assert(size(obj.loadVectorFunction([1;1;1]), 1) == 1, 'loadVectorFunction must be of size "1 x 1"');
                end
            else
                error('Not implemented yet!');
            end

            % check edof
            numberOfNodes = size(obj.meshObject.edof, 2);
            if strcmp(obj.loadGeometry, 'line')
                assert(numberOfNodes == obj.masterObject.shapeFunctionObject.order+1, ['edof must be of size "', num2str(obj.numberOfElements), ' x ', num2str(obj.masterObject.shapeFunctionObject.order+1), '"']);
            elseif strcmp(obj.loadGeometry, 'area')
                % TODO
            else
                error('Not implemented yet!');
            end
        end
    end
end
