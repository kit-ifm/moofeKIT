classdef dofClass < matlab.mixin.Copyable
    properties
        qR
        qN
        qN1
        vN
        vN1
        L
        J
        listContinuumObjects
        totalNumberOfDofs = 0;
        dirichletDof
        R
        K
        postDataObject
    end
    properties (Dependent = true)
        numberOfContinuumObjects
        solveDof
        doSolve
    end
    %% constructor
    methods
        function obj = dofClass
            obj.postDataObject = postDataClass();
        end
    end
    methods
        function initializeObjectsShapeFunctions(obj)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.initializeShapeFunctions;
            end
        end
        function initializeObjectsGlobalDofs(obj)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.initializeGlobalDofs(obj);
                if isa(continuumObject,'solidViscoClass')
                    continuumObject.initializeEpsilonViscoN1;
                end
            end
        end
        function initializeObjectsQV(obj)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.initializeQV;
            end
        end
        function initialize(obj,fieldNameCell)
            for index1 = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{index1};
                assert(ischar(fieldName));
                obj.(fieldName) = zeros(obj.totalNumberOfDofs,1);
            end
        end
        function updateGlobalField(obj,fieldNameCell)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateGlobalField(obj,fieldNameCell)
            end
        end
        function M = callMassMatrixElement(obj,setupObject)
            massDataFE = [];
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if continuumObject.callMassMatrix
                    if continuumObject.materialObject.rho > 0
                        dataFE = callElements(continuumObject, setupObject, 'massMatrix');
                    else
                        dataFE = continuumObject.storageFEObject.initializeDataFE('massMatrix');
                    end
                    massDataFE = [massDataFE, dataFE];
                end
            end
            % assembly
            M = sparse(vertcat(massDataFE(:).indexMeI),vertcat(massDataFE(:).indexMeJ),vertcat(massDataFE(:).Me),obj.totalNumberOfDofs,obj.totalNumberOfDofs);
        end
        function solveDof = get.solveDof(obj)
            solveDof = 1:obj.totalNumberOfDofs;
            if ~isempty(obj.dirichletDof)
                solveDof(obj.dirichletDof(:,1)) = [];
            end
        end
        function doSolve = get.doSolve(obj)
            doSolve = true(obj.totalNumberOfDofs,1);
            if ~isempty(obj.dirichletDof)
                doSolve(obj.dirichletDof) = false;
            end
        end
        function updateObjectsContinuumFieldPreNewtonLoop(obj,setupObject,fieldNameCell)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateContinuumFieldPreNewtonLoop(obj,fieldNameCell);
            end
        end
        function updateHistoryField(obj)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if ismethod(continuumObject,'updateHistoryField')
                    continuumObject.updateHistoryField(obj);
                end
            end
        end
        function updateTimeDependentFieldPreNewtonLoop(obj,fieldName,time) 
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateTimeDependentFieldPreNewtonLoop(obj,fieldName,time)
            end
        end
        function updateTimeDependentFieldNewtonLoop(obj,fieldName,time) 
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateTimeDependentFieldNewtonLoop(obj,fieldName,time)
            end
        end
        function updateObjectsContinuumFieldNewtonLoop(obj,setupObject,delta,fieldName)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateContinuumFieldNewtonLoop(obj,setupObject,delta,fieldName);
            end
        end
        function updateObjectsContinuumFieldPostNewtonLoop(obj,setupObject,fieldName)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                continuumObject.updateContinuumFieldPostNewtonLoop(obj,setupObject,fieldName);
            end
        end
        
        %%%%%%%%%%%%%%%%%
        function out = get.numberOfContinuumObjects(obj)
            out = numel(obj.listContinuumObjects);
        end
        function dirichletDof = get.dirichletDof(obj)
            dirichletDof = [];
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if ~isempty(continuumObject)
                    if isa(continuumObject,'dirichletClass')
                        dirichletDof = unique([dirichletDof; continuumObject.globalNodesDof]);
                    end
                end
            end
        end
    end
end
