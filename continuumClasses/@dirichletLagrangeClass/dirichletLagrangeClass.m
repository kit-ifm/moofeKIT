classdef dirichletLagrangeClass < baseFEClass
    properties (Constant)
        callMassMatrix = false;
        callElements = true;
    end
    properties
        storageFEObject
        masterObject                % object on which load is acting
        nodeList
        nodalDof = 1;               % 1,2,3         
        timeFunction = @(t) 0;      % time funktion of displacements
        time = 0;
        qR = [];                    % reference lagrange multipliers
        qN = [];                    % lagrange multipliers time step n
        qN1 = [];                   % lagrange multipliers time step n+1        
        vN = [];                    % velocity of lagrange multiplier time step n
        globalNodesDof
        dimension = 3;
    end

    %% constructor
    methods
        function obj = dirichletLagrangeClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
            obj.storageFEObject = storageFEClass();
        end
    end
    
    %% mandatory methods
    methods
        function initializeShapeFunctions(obj,dofObject)
            % nothing to do here
        end
        function initializeGlobalDofs(obj,dofObject)
            meshObject = obj.meshObject;
            masterMeshObject = obj.masterObject.meshObject;
            obj.nodeList = unique(obj.nodeList);
            numberOfDofs = numel(obj.nodeList)*numel(obj.nodalDof);
            meshObject.globalNodesDof = (1:numberOfDofs)' + dofObject.totalNumberOfDofs;
            globalFullEdof = zeros(numberOfDofs,2);
            temp = masterMeshObject.globalNodesDof(obj.nodeList,obj.nodalDof);
            globalFullEdof(:,1) = temp(:);
            globalFullEdof(:,2) = meshObject.globalNodesDof;
            meshObject.globalFullEdof = globalFullEdof;
            dofObject.totalNumberOfDofs = dofObject.totalNumberOfDofs + numberOfDofs;
        end
        function initializeQV(obj,dofObject)
            obj.qN = zeros(numel(obj.nodeList)*numel(obj.nodalDof),1);
            obj.qN1 = obj.qN;
        end
        function updateGlobalField(obj,dofObject,fieldNameCell)
            % nothing to do here
        end
        function updateContinuumFieldPreNewtonLoop(obj,dofObject,fieldNameCell)
            for index2 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index2};
                field = dofObject.(fieldName);
                assert(ischar(fieldName));
                if isprop(obj,fieldName)
                    obj.(fieldName) = field(obj.meshObject.globalNodesDof);
                end
            end
        end
        function updateTimeDependentFieldPreNewtonLoop(obj,dofObject,fieldName,time)
            % nothing to do here
        end
        function updateTimeDependentFieldNewtonLoop(obj,dofObject,fieldName,time)
            assert(ischar(fieldName));
            obj.time = time;
            if isprop(obj,fieldName)
                dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
            end
        end
        function updateContinuumFieldNewtonLoop(obj,dofObject,setupObject,delta,fieldName)
            field = dofObject.(fieldName);
            assert(ischar(fieldName));
            if isprop(obj,fieldName)
                obj.(fieldName) = field(obj.meshObject.globalNodesDof);
            end
        end
    end
    %% get und set methods
    methods
        function set.masterObject(obj,input)
            assert(isa(input,'solidClass')||isa(input,'solidThermoClass'),'masterObject must be compatibel to dirichlet class');
            obj.masterObject = input;
        end
        function set.timeFunction(obj,input)
            assert(isa(input,'function_handle'),'timeFunction must be of type function_handle');
            obj.timeFunction = input;
        end
        function set.nodeList(obj,input)
            checkIntegerMatrix(input);
            obj.nodeList = int32(input(:));
        end
        function set.nodalDof(obj,input)
            checkIntegerMatrix(input);
            obj.nodalDof = input;
        end
        function set.time(obj,t)
            assert(isnumeric(t),'time need to be of class double');
            obj.time = t; 
        end        
    end
end