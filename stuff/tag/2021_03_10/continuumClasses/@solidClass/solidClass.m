classdef solidClass < handle
    properties
        % mesh
        nodes
        edof                        % local edof
        fullEdof                    % local full edof
        orderShapeFunctions
        numberOfGausspoints
        dimension                       % 1,2,3
        shapeFunctions
        % material
        materialName
        materialData
        % solver
        globalFullEdof
        globalNodesDof
        qN
        qN1
        vN
        vN1
        zeroDirichletBC
        elementDisplacementType = 'displacement';
    end
    properties (Constant)
        flagNewton = struct('initializeGlobalDofs',true,...
                            'initializeQV',true,...
                            'initializeShapeFunctions',true,...
                            'updateTimeDependentField',false,...
                            'updateContinuumFieldPreNewtonLoop',true,...
                            'callElements',true,...
                            'updateContinuumFieldNewtonLoop',true);
    end
    properties (Dependent = true)
        qR
        massMatrix
    end
    methods
        %% constructor
        function obj = solidClass(dofsObject)
        end
        %% set and get methods
        function set.edof(obj,input)
            checkIntegerMatrix(input);
            obj.edof = int32(input);
        end
        function set.globalFullEdof(obj,input)
            checkIntegerMatrix(input);
            obj.globalFullEdof = int32(input);
        end
        function set.globalNodesDof(obj,input)
            checkIntegerMatrix(input);
            obj.globalNodesDof = int32(input);
        end
        %         function updateQV(obj,input)
        %             assert(ischar(input(1)) & isnumeric(input(2)));
        %             obj.(input(1)) = input(2);
        %         end
%         function updateField(obj,dofsObject,fieldNameCell,flagUpdate)
%             for index1 = 1:numel(fieldNameCell)
%                 fieldName = fieldNameCell{index1};
%                 field = dofsObject.(fieldName);
%                 assert(ischar(fieldName) & ischar(flagUpdate));
%                 for index2 = 1:numel(obj)
%                     if  isfield(obj(index2).flagNewton,flagUpdate)
%                         if obj(index2).flagNewton.(flagUpdate)
%                             obj(index2).(fieldName) = field(obj(index2).globalNodesDof);
%                         end
%                     end
%                 end
%             end
%         end
        function set.dimension(obj,input)
            checkIntegerMatrix(input);
            assert((input <= 3) & (input >= 1));
            obj.dimension = input;
        end
        function out = get.qR(obj)
            out = obj.nodes;
        end
        function plot(obj)
            for index = 1:numel(obj)
                patchData = struct('Vertices',[],'Faces',[]);
                %     element = '3Dbrick';
                SX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
                patchData.Vertices = obj(index).qN1;
                patchData.Faces = [ obj(index).edof(:,SX(1,:));
                    obj(index).edof(:,SX(2,:));
                    obj(index).edof(:,SX(3,:));
                    obj(index).edof(:,SX(4,:));
                    obj(index).edof(:,SX(5,:));
                    obj(index).edof(:,SX(6,:))];
                patchData.FaceColor = [0 1 0];
                patchData.FaceAlpha = 0.5;
                patch(patchData);
                hold on
            end
        end
        function out = get.massMatrix(obj)
            out = massMatrixElement(obj);
        end
                %% helpful methods
        %        function checkIntegerMatrix2(obj,input)
        %            assert(isnumeric(input) & (sum(sum(mod(input,1))) == 0),'input need to be of type integer even if it is classified as float');
        %        end
        %% solver initialization
        out = assignGlobalDofs(obj,input)
        initializeQV(obj,input)
        assignShapefunctions(obj)
    end
end