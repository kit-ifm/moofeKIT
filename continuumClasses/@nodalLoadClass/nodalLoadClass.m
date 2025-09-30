classdef nodalLoadClass < baseFEClass
    %NODALLOADCLASS Nodal loads.

    %% mandatory properties
    properties (Constant)
        callMassMatrix = false;
        callElements = false;
    end

    %% further properties
    properties
        masterObject
        loadVector
        loadType = 'deadLoad'; % 'deadLoad' or 'followerLoad'
        loadPhysics = 'mechanical'; % 'mechanical', 'thermal', 'electrical', tbc
        timeFunction = @(t) 1;
        nodeList
    end
    properties (SetAccess = private)
        elementData = struct();
    end

    %% constructor
    methods
        function obj = nodalLoadClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end

    %% mandatory methods
    methods
        function initializeShapeFunctions(obj, dofObject)
            % nothing to do here
        end
        function initializeGlobalDofs(obj, dofObject)
            % nothing to do here
        end
        function initializeQV(obj, dofObject)
            % nothing to do here
            verificationOfInputData(obj)
        end
        function updateGlobalField(obj, dofObject, fieldNameCell)
            % nothing to do here
        end
        function updateContinuumFieldPreNewtonLoop(obj, setupObject, fieldNameCell)
            % nothing to do here
        end
        function updateTimeDependentFieldPreNewtonLoop(obj, dofObject, fieldName, time)
            % nothing to do here
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

    %% set methods
    methods

        function set.masterObject(obj, input)
            assert(isa(input, 'solidSuperClass'), 'masterObject of nodalForce object must be of type solidSuperClass!')
            obj.masterObject = input;
        end

        function set.loadVector(obj, value)
            assert(size(value, 1) == 1 || size(value, 2) == 1, 'loadVector must be a vector!');
            if size(value, 1) == 1
                value = value.';
            end
            obj.loadVector = value;
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

        function set.nodeList(obj, input)
            checkIntegerMatrix(input);
            obj.nodeList = int32(input(:));
        end
    end

    %% verification of input data
    methods
        function verificationOfInputData(obj)
            % load vector
            if strcmp(obj.loadPhysics, 'mechanical')
                if isa(obj.masterObject,'beamClass')
                    displacementDofsPerNode = size(obj.masterObject.qR, 2);
                else
                    displacementDofsPerNode = size(obj.masterObject.qR, 2) - obj.masterObject.additionalFields;
                end
                assert(size(obj.loadVector, 1) == displacementDofsPerNode, ['loadVector must be of size "', num2str(displacementDofsPerNode), ' x 1"']);
            elseif any(strcmp(obj.loadPhysics, {'thermal', 'electrical'}))
                assert(size(obj.loadVector, 1) == 1, 'loadVector must be of size "1 x 1"');
            else
                error('Not implemented yet!');
            end
        end
    end

    %% compute dataFE
    methods
        function R = computeNodalForceContribution(obj, dofObject, setupObject)
            assert(~strcmp(obj.loadType, 'followerLoad'), 'loadType="followerLoad" not implemented yet!');
            R = zeros(dofObject.totalNumberOfDofs, 1);
            globalNodesDof = obj.masterObject.meshObject.globalNodesDof;

            if isa(obj.masterObject,'beamClass')
                displacementDofsPerNode = size(obj.masterObject.qR, 2);
            else
                displacementDofsPerNode = size(obj.masterObject.qR, 2) - obj.masterObject.additionalFields;
            end

            if strcmp(obj.loadPhysics, 'mechanical')
                globalNodesDof = globalNodesDof(:, 1:displacementDofsPerNode);
            elseif strcmp(obj.loadPhysics, 'thermal')
                if isa(obj.masterObject, 'solidThermoClass')
                    globalNodesDof = globalNodesDof(:, displacementDofsPerNode+1);
                elseif isa(obj.masterObject, 'solidElectroThermoClass')
                    globalNodesDof = globalNodesDof(:, displacementDofsPerNode+2);
                else
                    error('Not implemented yet!');
                end
            elseif strcmp(obj.loadPhysics, 'electrical')
                if isa(obj.masterObject, 'solidElectroClass') || isa(obj.masterObject, 'solidElectroThermoClass')
                    globalNodesDof = globalNodesDof(:, displacementDofsPerNode+1);
                else
                    error('Not implemented yet!');
                end
            end

            timeEvaluationPoint = setupObject.time;
            if any(strcmp(setupObject.integrator, {'Midpoint', 'DiscreteGradient'}))
                timeEvaluationPoint = timeEvaluationPoint - 1 / 2 * setupObject.timeStepSize;
            end
            for i=1:length(obj.nodeList)
                R(globalNodesDof(obj.nodeList(i), :)) = obj.loadVector * obj.timeFunction(timeEvaluationPoint);
            end
        end

        function computeExternalEnergyContribution(obj, dofObject, setupObject)
            numberOfNodes = length(obj.nodeList);
            if isa(obj.masterObject,'beamClass')
                displacementDofsPerNode = size(obj.masterObject.qR, 2);
            else
                displacementDofsPerNode = size(obj.masterObject.qR, 2) - obj.masterObject.additionalFields;
            end
            fieldsToConsider = 1:displacementDofsPerNode;
            if strcmp(obj.loadPhysics, 'thermal')
                if isa(obj.masterObject, 'solidThermoClass')
                    fieldsToConsider = displacementDofsPerNode + 1;
                elseif isa(obj.masterObject, 'solidElectroThermoClass')
                    fieldsToConsider = displacementDofsPerNode + 2;
                else
                    error('Not implemented yet!');
                end
            elseif strcmp(obj.loadPhysics, 'electrical')
                if isa(obj.masterObject, 'solidElectroClass') || isa(obj.masterObject, 'solidElectroThermoClass')
                    fieldsToConsider = displacementDofsPerNode + 1;
                else
                    error('Not implemented yet!');
                end
            end
            uN1 = (obj.masterObject.qN1(obj.nodeList, fieldsToConsider) - obj.masterObject.qR(obj.nodeList, fieldsToConsider)).'; % displacement of nodes
            F = repmat(obj.loadVector, numberOfNodes, 1)* obj.timeFunction(setupObject.time); % force on nodes
            obj.elementData(setupObject.timeStep+1).externalEnergy = uN1(:).' * F;
            if setupObject.timeStep ==1 %% initialise for t=0
                obj.elementData(1).externalEnergy = 0;
            end
        end
    end
end
