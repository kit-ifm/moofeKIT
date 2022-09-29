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
        listOfAllContinuumObjects
        numberOfContinuumObjects
        callMassMatrixElement
        solveDof
    end
    methods
        % constructor
        function obj = dofClass
            obj.collectSolverObjects;
        end
    end
    methods        
        function dirichletDof = get.dirichletDof(obj)
            dirichletDof = [];
            for index = 1:obj.numberOfContinuumObjects
                continuumObj = evalin('base',obj.listContinuumObjects{index});
                if ~isempty(continuumObj)
                    if isprop(continuumObj(1),'globalDirichletList')
                        dirichletDof = unique(vertcat(continuumObj.globalDirichletList));
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
        function out = get.listOfAllContinuumObjects(obj)
            listing = dir('../../continuumClasses');
            out = {};
            for index = 3:numel(listing)
                out{index-2} = strcat(listing(index).name(2:end-5),'Objects');
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
            for index = 1:obj.numberOfContinuumObjects
                continuumObj = evalin('base',obj.listContinuumObjects{index});
                for index1 = 1:numel(fieldNameCell)
                    fieldName = fieldNameCell{index1};
                    assert(ischar(fieldName));
                    for index2 = 1:numel(continuumObj)
                        if isprop(continuumObj(index2),fieldName)
                            obj.(fieldName)(continuumObj(index2).globalNodesDof) = continuumObj(index2).(fieldName);
                        end
                    end
                end
            end
        end
        function updateTimeDependentField(obj,fieldNameCell,time,flagUpdate) 
            for index = 1:obj.numberOfContinuumObjects
                continuumObj = evalin('base',obj.listContinuumObjects{index});                
                fieldName = fieldNameCell{1};
                assert(ischar(fieldName));
                for index1 = 1:numel(continuumObj)
                    if isprop(continuumObj(index1),fieldName)
                        if continuumObj(index1).flagNewton.(flagUpdate)
                            continuumObj(index1).time = time;
                            obj.(fieldName)(continuumObj(index1).globalDirichletList) = continuumObj(index1).(fieldName);
                        end
                    end
                end
            end
        end        
        function updateContinuumField(obj,fieldNameCell,flagUpdate)
            for index = 1:obj.numberOfContinuumObjects
                continuumObj = evalin('base',obj.listContinuumObjects{index});
                for index1 = 1:numel(fieldNameCell)
                    fieldName = fieldNameCell{index1};
                    field = obj.(fieldName);
                    assert(ischar(fieldName) & ischar(flagUpdate));
                    for index2 = 1:numel(continuumObj)
                        if isprop(continuumObj(index2),fieldName) && isfield(continuumObj(index2).flagNewton,flagUpdate)
                            if continuumObj(index2).flagNewton.(flagUpdate)
                                continuumObj(index2).(fieldName) = field(continuumObj(index2).globalNodesDof);
                            end
                        end
                    end
                end
            end
        end
        function out = get.callMassMatrixElement(obj)
            massData = [];
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = evalin('base',obj.listContinuumObjects{index1});
                for index2 = 1:numel(continuumObject)
                    if isprop(continuumObject(index2),'massMatrix')
                        massData = [massData, continuumObject(index2).massMatrix];
                    end
                end
            end
            % assembly
            out = sparse(vertcat(massData(:).indexMi),vertcat(massData(:).indexMj),vertcat(massData(:).MeVector),obj.totalNumberOfDofs,obj.totalNumberOfDofs);
        end
        function initializeGlobalDofs(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = evalin('base',obj.listContinuumObjects{index1});
                for index2 = 1:numel(continuumObject)
                    if isfield(continuumObject(index2).flagNewton,flag)
                        if continuumObject(index2).flagNewton.(flag)
                            numberOfDofPerNode = size(continuumObject(index2).nodes,2);
                            edof = continuumObject(index2).edof;
                            numberOfNodes = size(continuumObject(index2).nodes,1);
                            numberOfElements = size(edof,1);
                            globalFullEdof = zeros(numberOfElements,size(edof,2)*numberOfDofPerNode);
                            for j = 1:numberOfDofPerNode
                                globalFullEdof(:,j:numberOfDofPerNode:end) = obj.totalNumberOfDofs + edof*numberOfDofPerNode-(numberOfDofPerNode-j);
                            end
                            globalNodesDof = zeros(numberOfNodes,numberOfDofPerNode);
                            globalNodesDof(:,1:numberOfDofPerNode) = obj.totalNumberOfDofs + kron(ones(numberOfNodes,1),1:numberOfDofPerNode)+kron((0:numberOfDofPerNode:numberOfDofPerNode*numberOfNodes-1)',ones(1,numberOfDofPerNode));
                            continuumObject(index2).globalFullEdof = globalFullEdof;
                            continuumObject(index2).globalNodesDof = globalNodesDof;
                            obj.totalNumberOfDofs = obj.totalNumberOfDofs + numberOfNodes*numberOfDofPerNode;
                        end
                    end
                end
            end
        end
        function initializeQV(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = evalin('base',obj.listContinuumObjects{index1});
                for index2 = 1:numel(continuumObject)
                    if isfield(continuumObject(index2).flagNewton,flag)
                        if continuumObject(index2).flagNewton.(flag)
                            if isempty(continuumObject(index2).qN)
                                continuumObject(index2).qN = continuumObject(index2).qR;
                            end
                            if isempty(continuumObject(index2).qN1)
                                continuumObject(index2).qN1 = continuumObject(index2).qR;
                            end
                            if isempty(continuumObject(index2).vN)
                                continuumObject(index2).vN = zeros(size(continuumObject(index2).nodes,1), size(continuumObject(index2).nodes,2));
                            end
                        end
                    end
                end
            end
        end
        function initializeShapeFunctions(obj,flag)
            for index1 = 1:obj.numberOfContinuumObjects
                continuumObject = evalin('base',obj.listContinuumObjects{index1});
                for index2 = 1:numel(continuumObject)
                    if isfield(continuumObject(index2).flagNewton,flag)
                        if continuumObject(index2).flagNewton.(flag)
                            continuumObject(index2).shapeFunctions = lagrangeShapeFunctions(size(continuumObject(index2).edof,2), continuumObject(index2).numberOfGausspoints, continuumObject(index2).dimension);
                        end
                    end
                end
            end
        end
        function obj = collectSolverObjects(obj)
            workSpace = evalin('base','whos');
            namesOfWorkSpace = {workSpace.name}.';
            classesOfWorkSpace = {workSpace.class}.';
            indexList = 0;
            for index1 = 1:numel(obj.listOfAllContinuumObjects)
                objectsName = obj.listOfAllContinuumObjects{index1}(1:end-7);
                objectsClassName = strcat(objectsName,'Class');
                listClassObjectsWorkSpace = namesOfWorkSpace(strcmp(classesOfWorkSpace,objectsClassName));               
                continuumObjectsCell = {};
                if numel(listClassObjectsWorkSpace) > 0
                    index2 = 1;
                    for index3 = 1:numel(listClassObjectsWorkSpace)
                        continuumObjectsWorkSpace = evalin('base',listClassObjectsWorkSpace{index3});
                        if isa(continuumObjectsWorkSpace,objectsClassName)
                            for index4 = 1:numel(continuumObjectsWorkSpace)
                                continuumObjectsCell{index2} = continuumObjectsWorkSpace(index4);
                                index2 = index2 + 1;
                            end
                        else
                            error('Class of object not found!')
                        end
                    end
                end
                if numel(continuumObjectsCell) > 0
                    continuumObjects(1,numel(continuumObjectsCell)) = eval(objectsClassName);
                    for index3 = 1:numel(continuumObjectsCell)
                        continuumObjects(index3) = continuumObjectsCell{index3};
                    end
                    assignin('base',strcat(objectsName,'Objects'), continuumObjects)
                    clear continuumObjects
                    indexList = indexList + 1;
                    obj.listContinuumObjects{indexList} = strcat(objectsName,'Objects');
                end
            end
        end        
    end
end
