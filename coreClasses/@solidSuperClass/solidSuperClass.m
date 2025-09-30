classdef solidSuperClass < baseFEClass
% SOLIDSUPERCLASS Super class inherited from baseFEClass for all finite 
% element classes for continuum solid objects which defines necessary 
% methods for all child classes.
% 
% INHERIT
% Inherit a child class via classdef childClass < SOLIDSUPERCLASS. 
% 
% GENERAL
% A superclass class serves as a basis for a group of related subclasses. 
% A superclass can define abstract properties and methods that subclasses 
% implement. Each subclass can implement the concrete properties and 
% methods in a way that supports their specific requirements.
% 
% REFERENCE(S)
% https://de.mathworks.com/help/matlab/matlab_oop/calling-superclass-methods-on-subclass-objects.html 
% https://de.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html
% 
% SEE ALSO
% baseFEClass, 
% solidClass, 
% solidThermoClass, 
% solidElectroThermoClass 
% 
% CREATOR(S) 
% Marlon Franke

    properties (Abstract = true)
        additionalFields
        dofsPerAdditionalField
    end
    properties (Constant)
        callMassMatrix = true;
        callElements = true;
    end
    properties
        dimension
        elementDisplacementType = 'displacement';   % mixedSC
        elementNameAdditionalSpecification = '';    % for mixedSC: NonCascadeCGc, InvariantCGc; attention: empty '' means cascade version of CGc
        qR
        elementData = struct();
        momenta = struct('L',[],'J',[]);
        %flagNumericalTangent = [false;...   % compute element with numerical tangent
        %                        false];     % show difference of numerical and analytical tangent
        % post
        plotData = struct('handleMesh',[],'handlePatch',[],'surfElements',[],'surfElementsKey',[],'N',[],'GeoData',[],'nodalColorData',[]);
        L = [];
        % solver
        vN
        vN1
        qN
        qN1
    end
    properties (Dependent = true)
        J
    end
    properties (SetAccess=private)
        elementGeometryType
    end
    properties
        % load objects
        shapeFunctionObject% = lagrangeShapeFunctionClass();
        storageFEObject% = storageFEClass();
        materialObject% = storageFEClass();
        mapVoigtObject% = mapVoigtClass();
        mixedFEObject% = mixedFEClass();
        numericalTangentObject
        artificialNeuralNetworkObject
    end
    methods
        function obj = solidSuperClass()
        % constructor
            obj.meshObject = meshClass();
            obj.storageFEObject = storageFEClass();
            obj.shapeFunctionObject = shapeFunctionClass();
            obj.mapVoigtObject = mapVoigtClass();
            obj.mixedFEObject = mixedFEClass();
            obj.numericalTangentObject = numericalTangentClass();
            obj.artificialNeuralNetworkObject = artificialNeuralNetworkClass();
        end
    end
    methods(Access = protected)
        % override copyElement method:
        function copyObject = copyElement(obj)
            % make a shallow copy of all properties
            copyObject = copyElement@matlab.mixin.Copyable(obj);
            % insert all subObjects which should provide a deep copy behaviour
            copyObject.meshObject = copy(obj.meshObject);
            copyObject.storageFEObject = copy(obj.storageFEObject);
            copyObject.shapeFunctionObject = copy(obj.shapeFunctionObject);
            copyObject.mapVoigtObject = copy(obj.mapVoigtObject);
            copyObject.mixedFEObject = copy(obj.mixedFEObject);
            copyObject.numericalTangentObject = copy(obj.numericalTangentObject);
        end
    end
    methods
        function out = get.J(obj)                          % Get method for the angular momentum.
            if isa(obj,"beamClass")
                out = obj.getAngularMomentum;
            else
                out = 0;
                if size(obj.qN1,2) == 0 || size(obj.L,2) == 0
                    if obj.dimension == 2
                        out = 0;
                    elseif obj.dimension == 3
                        out = [0 0 0]';
                    end
                else
                    if obj.dimension == 2
                        out = sum(obj.qN1(:,1).*obj.L(:,2) - obj.qN1(:,2).*obj.L(:,1));
                    elseif obj.dimension == 3
                        out = [ sum(obj.qN1(:,2).*obj.L(:,3) - obj.qN1(:,3).*obj.L(:,2)),...
                                sum(obj.qN1(:,3).*obj.L(:,1) - obj.qN1(:,1).*obj.L(:,3)),...
                                sum(obj.qN1(:,1).*obj.L(:,2) - obj.qN1(:,2).*obj.L(:,1))];
                    end
                end
            end
        end
        function rotate(obj,rotatingAxis,alpha)
            assert(obj.dimension==3,'the rotating operation is only implemented for 3D')
            assert(rotatingAxis==1 || rotatingAxis==2 || rotatingAxis==3,'the chosen rotating axis must be 1, 2 or 3')
            switch rotatingAxis
                case 1
                    R = [   1 0 0;...
                            0 cos(alpha) -sin(alpha);...
                            0 sin(alpha) cos(alpha)];
                case 2
                    R = [   cos(alpha)  0 -sin(alpha);...
                            0          1           0;...
                            sin(alpha)  0 cos(alpha)];
                case 3
                    R = [   cos(alpha) -sin(alpha) 0;...
                            sin(alpha) cos(alpha) 0;...
                            0 0 1];
            end
            obj.meshObject.nodes(:,1:3) = (R*(obj.meshObject.nodes(:,1:3))')';
        end
        function reflection(obj,reflectionAxis)
            assert(obj.dimension==3,'the reflection operation is only implemented for 3D')
            assert(reflectionAxis==1 || reflectionAxis==2 || reflectionAxis==3,'the chosen reflection axis must be 1, 2 or 3')
            switch reflectionAxis
                case 1
                    R = [  -1 0 0;...
                            0 1 0;...
                            0 0 1];
                case 2
                    R = [   1 0 0;...
                            0 -1 0;...
                            0 0 1];
                case 3
                    R = [   1 0 0;...
                            0 1 0;...
                            0 0 -1];
            end
            obj.meshObject.nodes(:,1:3) = (R*(obj.meshObject.nodes(:,1:3))')';
            obj.qR(:,1:3) = (R*(obj.qR(:,1:3))')';
            obj.qN(:,1:3) = (R*(obj.qN(:,1:3))')';
            obj.qN1(:,1:3) = (R*(obj.qN1(:,1:3))')';
        end
    end
    %% mandatory methods
    methods 
        plot(obj,setupObject)
        function initializeShapeFunctions(obj)
        % method                    initializeShapeFunctions
        % mandatory implementation  inherited abstract method from baseFEClass
            numberOfNodes = size(obj.meshObject.edof,2);

            % determine elementGeometryType
            obj.elementGeometryType = determineElementGeometryType(obj.dimension, numberOfNodes);

            % compute shape functions
            obj.shapeFunctionObject.computeShapeFunction(obj.dimension, numberOfNodes, obj.elementGeometryType);
            if isprop(obj,'selectiveReducedShapeFunctionObject')
                if isempty(obj.selectiveReducedShapeFunctionObject.order) 
                    obj.selectiveReducedShapeFunctionObject.order = obj.shapeFunctionObject.order;
                end
                if isempty(obj.selectiveReducedShapeFunctionObject.numberOfGausspoints)
                    obj.selectiveReducedShapeFunctionObject.numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                end

                obj.selectiveReducedShapeFunctionObject.computeShapeFunction(obj.dimension, numberOfNodes, obj.elementGeometryType);
                
            end
        end
        function initializeGlobalDofs(obj,dofObject)
        % method                    initializeGlobalDofs
        % mandatory implementation  inherited abstract method from baseFEClass
            %numberOfDofPerNode = size(obj.meshObject.nodes,2);
            numberOfDofPerNode = size(obj.qR,2);
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
            if ~contains(obj.elementDisplacementType,{'displacement', 'thermo'})
                obj.mixedFEObject.initializeMixedElements(dofObject,obj);
                obj.meshObject.globalFullEdof = [obj.meshObject.globalFullEdof, obj.mixedFEObject.globalEdof];
            end
        end
        function initializeQR(obj)
            % method initializeQR
            obj.qR = obj.meshObject.nodes;
        end
        function initializeQV(obj)
        % method                    initializeQV
        % mandatory implementation  inherited abstract method from baseFEClass
            obj.initializeQR;
            if isempty(obj.qN)
                obj.qN = obj.qR;
            end
            if isempty(obj.qN1)
                obj.qN1 = obj.qN;
            end
            if isempty(obj.vN)
                obj.vN = zeros(size(obj.qR));
            end
            if isempty(obj.vN1)
                obj.vN1 = zeros(size(obj.meshObject.nodes,1), size(obj.meshObject.nodes,2));
            end
        end
        function updateGlobalField(obj,dofObject,fieldNameCell)
        % method                    updateGlobalField
        % mandatory implementation  inherited abstract method from baseFEClass
            for index2 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index2};
                assert(ischar(fieldName));
                if isprop(obj,fieldName)
                    dofObject.(fieldName)(obj.meshObject.globalNodesDof) = obj.(fieldName);
                end
            end
            updateGlobalField(obj.mixedFEObject,obj,dofObject);

        end
        function updateContinuumFieldPreNewtonLoop(obj,dofObject,fieldNameCell)
        % method                    updateContinuumFieldPreNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
            for index2 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index2};
                field = dofObject.(fieldName);
                assert(ischar(fieldName));
                if isprop(obj,fieldName)
                    obj.(fieldName) = field(obj.meshObject.globalNodesDof);
                end
            end
            % mixed elements update
            if ~strcmpi(obj.elementDisplacementType,'displacement')
                if ~isempty(obj.mixedFEObject)
                    obj.mixedFEObject.qN = obj.mixedFEObject.qN1;
                end
            end
        end
        function updateTimeDependentFieldPreNewtonLoop(obj,dofObject,fieldName,time)
        % method                    updateTimeDependentFieldPreNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
        % nothing to do for solids
        end
        function updateTimeDependentFieldNewtonLoop(obj,dofObject,fieldName,time)
        % method                    updateTimeDependentFieldNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
        % nothing to do for solids
        end
        function updateContinuumFieldNewtonLoop(obj,dofObject,setupObject,delta,fieldName)
        % method                    updateContinuumFieldNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
            field = dofObject.(fieldName);
            assert(ischar(fieldName));
            if isprop(obj,fieldName)
                obj.(fieldName) = field(obj.meshObject.globalNodesDof);
                % mixed elements:
                if isa(obj,'solidSuperClass') && ~strcmpi(obj.elementDisplacementType,'displacement')
                    if ~isempty(obj.mixedFEObject)
                        obj.mixedFEObject.updateMixedElements(obj,dofObject,setupObject,delta,fieldName,field)
                    end
                end
            end
        end
        function updateContinuumFieldPostNewtonLoop(obj,dofObject,setupObject,fieldName)
        % method                    updateContinuumFieldPostNewtonLoop
        % mandatory implementation  inherited abstract method from baseFEClass
            field = dofObject.(fieldName);
            assert(ischar(fieldName));
            if isprop(obj,fieldName)
                obj.(fieldName) = field(obj.meshObject.globalNodesDof);
            end
        end
    end
    %% get and set methods
    methods
        function set.qR(obj, value)
            % set method qR
            assert(ismatrix(value))
            obj.qR = value;
        end

        function set.elementNameAdditionalSpecification(obj, value)
            % set method elementNameAdditionalSpecification
            assert(ischar(value), "Property 'elementNameAdditionalSpecification' must be of type string!");
            obj.elementNameAdditionalSpecification = value;
        end
    end
end