classdef neumannClass < baseFEClass
    properties (Constant)
        callMassMatrix = false;
        callElements = true;
    end
    properties
        storageFEObject
        shapeFunctionObject% = lagrangeShapeFunctionClass();
        masterObject
        %         edof
        forceVector
        typeOfLoad = 'deadLoad'       % deadLoad or followerLoad
        time = 0;
        timeFunction = @(t) 1;
        %TODO
        projectionType = 'none'       % which length is used for calculation of area (not needed in 1D) --> 'none': real length, 'x': length in xdirection, 'y': length in ydirection
        globalEdof
        shapeFunctions
        masterNodes
        nodalForce
        dimension = 3;
        field = 'mechanical'           % mechanical, thermal, electrical, tbc
        omega_0
        energy

        ePot = struct('externalEnergy',0);
    end
    properties (SetAccess=private)
        elementGeometryType
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
        %
        plot(obj, setupObject)
        %
        function initializeShapeFunctions(obj, dofObject)
            dimensionMasterObject = obj.masterObject.dimension;
            if dimensionMasterObject == 1
                % no shape functions in 1D as only nodal loads
            elseif dimensionMasterObject == 2 || dimensionMasterObject == 3
                numberOfNodes = size(obj.meshObject.edof,2);
                obj.elementGeometryType = determineElementGeometryType(dimensionMasterObject-1, numberOfNodes);
                obj.shapeFunctionObject.computeShapeFunction(dimensionMasterObject-1, numberOfNodes, obj.elementGeometryType);                
            end
        end
        function initializeGlobalDofs(obj, dofObject)
            edof = obj.meshObject.edof;

            if strcmpi(obj.field, 'mechanical')
                globalEdof = zeros(obj.numberOfElements, size(edof, 2)*obj.dimension);
                for index2 = 1:obj.dimension
                    temp = zeros(obj.numberOfElements, size(edof, 2));
                    temp(:) = obj.masterObject.meshObject.globalNodesDof(edof, index2);
                    globalEdof(:, index2:obj.dimension:end) = temp;
                end
            elseif strcmpi(obj.field, 'thermal')
                if isa(obj.masterObject, 'solidThermoClass')
                    globalEdof = zeros(obj.numberOfElements, size(edof, 2));
                    globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, obj.dimension+1);
                elseif isa(obj.masterObject, 'solidElectroThermoClass')
                    globalEdof = zeros(obj.numberOfElements, size(edof, 2));
                    globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, obj.dimension+2);
                end
            elseif strcmpi(obj.field, 'electrical')
                globalEdof = zeros(obj.numberOfElements, size(edof, 2));
                globalEdof(:) = obj.masterObject.meshObject.globalNodesDof(edof, obj.dimension+1);
            end
            obj.globalEdof = globalEdof;
            nodes = unique(edof');
            obj.masterNodes = nodes;
            obj.nodalForce = zeros(size(nodes, 1), obj.masterObject.dimension);
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
            obj.time = time;
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
        function set.masterObject(obj, input)
            assert(isa(input, 'solidSuperClass'), 'masterObject of Neumann Object must be of type solidSuperClass!')
            obj.masterObject = input;
        end
        function set.projectionType(obj, projection)
            if ~strcmp(projection, 'none') && ~strcmp(projection, 'x') && ~strcmp(projection, 'y')
                error('projection must be "none", "x" or "y"')
            else
                obj.projectionType = projection;
            end
        end
        function set.time(obj, time)
            obj.time = time;
        end
        function set.typeOfLoad(obj, type)
            if ~strcmp(type, 'deadLoad') && ~strcmp(type, 'followerLoad')
                error('loadtype must be either "deadLoad" or "followerLoad"')
            elseif strcmp(type, 'followerLoad')
                error('followerLoad loads currently not implemented')
            else
                obj.typeOfLoad = type;
            end
        end
        function set.forceVector(obj, forceVector)
            assert(size(forceVector, 1) == obj.masterObject.dimension, 'Force-vector must be of size "dim x 1"') %#ok<MCSUP>
            assert(size(forceVector, 2) == 1, 'Force-vector must be of size "dim x 1"')
            obj.forceVector = forceVector;
        end
        function setTime(obj, t)
            obj.time = t;
        end
        function out = get.numberOfElements(obj)
            out = size(obj.meshObject.edof, 1);
        end
        function updateNodalForce(obj)
            R = full(sparse(vertcat(obj.storageFEObject.dataFE(:).indexReI), 1, vertcat(obj.storageFEObject.dataFE(:).Re), max(vertcat(obj.storageFEObject.dataFE(:).indexReI)), 1));
            obj.nodalForce(:, :) = -R(obj.masterObject.meshObject.globalNodesDof(obj.masterNodes, :));
        end
    end
end
