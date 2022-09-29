classdef dofClass < handle
    properties
        qN
        qN1
        vN
        vN1
        listContinuumObjects
        totalNumberOfDofs = 0;
        dirichletDof
    end
    properties (Dependent = true)
        numberOfContinuumObjects
        callMassMatrixElement
        solveDof
    end
    methods
        % constructor
        function obj = dofClass
        end
    end
    methods        
        function dirichletDof = get.dirichletDof(obj)
            dirichletDof = [];
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if ~isempty(continuumObject)
                    if isprop(continuumObject,'globalDirichletList')
                        dirichletDof = unique([dirichletDof; continuumObject.globalDirichletList]);
                    end
                end
            end
        end
        function solveDof = get.solveDof(obj)
            solveDof = 1:obj.totalNumberOfDofs;
            if ~isempty(obj.dirichletDof)
                solveDof(obj.dirichletDof(:,1)) = [];
            end
        end    
        function out = get.numberOfContinuumObjects(obj)
            out = numel(obj.listContinuumObjects);
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
                for index2 = 1:numel(fieldNameCell)
                    fieldName = fieldNameCell{index2};
                    assert(ischar(fieldName));
                    if isprop(continuumObject,fieldName)
                        obj.(fieldName)(continuumObject.globalNodesDof) = continuumObject.(fieldName);
                    end
                end
            end
        end
        function updateTimeDependentField(obj,fieldNameCell,time,flagUpdate) 
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                fieldName = fieldNameCell{1};
                assert(ischar(fieldName));
                if isprop(continuumObject,fieldName)
                    if continuumObject.flagNewton.(flagUpdate)
                        continuumObject.time = time;
                        obj.(fieldName)(continuumObject.globalDirichletList) = continuumObject.(fieldName);
                    end
                end
            end
        end        
        function updateContinuumField(obj,fieldNameCell,flagUpdate)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                for index2 = 1:numel(fieldNameCell)
                    fieldName = fieldNameCell{index2};
                    field = obj.(fieldName);
                    assert(ischar(fieldName) & ischar(flagUpdate));
                    if isprop(continuumObject,fieldName) && isfield(continuumObject.flagNewton,flagUpdate)
                        if continuumObject.flagNewton.(flagUpdate)
                            continuumObject.(fieldName) = field(continuumObject.globalNodesDof);
                        end
                    end
                end
            end
        end
        function out = get.callMassMatrixElement(obj)
            massData = [];
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if isprop(continuumObject,'massMatrix')
                    massData = [massData, continuumObject.massMatrix];
                end
            end
            % assembly
            out = sparse(vertcat(massData(:).indexMi),vertcat(massData(:).indexMj),vertcat(massData(:).MeVector),obj.totalNumberOfDofs,obj.totalNumberOfDofs);
        end
        function initializeGlobalDofs(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if isfield(continuumObject.flagNewton,flag)
                    if continuumObject.flagNewton.(flag)
                        numberOfDofPerNode = size(continuumObject.nodes,2);
                        edof = continuumObject.edof;
                        numberOfNodes = size(continuumObject.nodes,1);
                        numberOfElements = size(edof,1);
                        globalFullEdof = zeros(numberOfElements,size(edof,2)*numberOfDofPerNode);
                        for j = 1:numberOfDofPerNode
                            globalFullEdof(:,j:numberOfDofPerNode:end) = obj.totalNumberOfDofs + edof*numberOfDofPerNode-(numberOfDofPerNode-j);
                        end
                        globalNodesDof = zeros(numberOfNodes,numberOfDofPerNode);
                        globalNodesDof(:,1:numberOfDofPerNode) = obj.totalNumberOfDofs + kron(ones(numberOfNodes,1),1:numberOfDofPerNode)+kron((0:numberOfDofPerNode:numberOfDofPerNode*numberOfNodes-1)',ones(1,numberOfDofPerNode));
                        continuumObject.globalFullEdof = globalFullEdof;
                        continuumObject.globalNodesDof = globalNodesDof;
                        obj.totalNumberOfDofs = obj.totalNumberOfDofs + numberOfNodes*numberOfDofPerNode;
                    end
                end
            end
        end
        function initializeQV(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if isfield(continuumObject.flagNewton,flag)
                    if continuumObject.flagNewton.(flag)
                        if isempty(continuumObject.qN)
                            continuumObject.qN = continuumObject.qR;
                        end
                        if isempty(continuumObject.qN1)
                            continuumObject.qN1 = continuumObject.qR;
                        end
                        if isempty(continuumObject.vN)
                            continuumObject.vN = zeros(size(continuumObject.nodes,1), size(continuumObject.nodes,2));
                        end
                    end
                end
            end
        end
        function initializeShapeFunctions(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = obj.listContinuumObjects{index1};
                if isfield(continuumObject.flagNewton,flag)
                    if continuumObject.flagNewton.(flag)
                        continuumObject.shapeFunctions = lagrangeShapeFunctions(size(continuumObject.edof,2), continuumObject.numberOfGausspoints, continuumObject.dimension);
                    end
                end
            end
        end
    end
end
