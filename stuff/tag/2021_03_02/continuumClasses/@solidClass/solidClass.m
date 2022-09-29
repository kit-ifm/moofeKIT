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
        globalNodesDOF
        QN
        QN1
        VN
        zeroDirichletBC
        elementDisplacementType = 'displacement';        
    end
    properties (Dependent = true)
        QR
    end    
    methods
        %% constructor
        function obj = solidClass
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
        function set.globalNodesDOF(obj,input)
            checkIntegerMatrix(input);
            obj.globalNodesDOF = int32(input);
        end
        function updateQV(obj,input)
            assert(ischar(input(1)) & isnumeric(input(2)));
            obj.(input(1)) = input(2);
        end
        function updateField(obj,fieldName,field)
            assert(ischar(fieldName) & isnumeric(field));
            for index = 1:numel(obj)
                obj(index).(fieldName) = field(obj(index).globalNodesDOF);
            end
        end
        function set.dimension(obj,input)
            checkIntegerMatrix(input);
            assert((input <= 3) & (input >= 1));
            obj.dimension = input;
        end
        function out = get.QR(obj)
            out = obj.nodes;
        end
        function plot(obj)
            for index = 1:numel(obj)
                patchData = struct('Vertices',[],'Faces',[]);
                %     element = '3Dbrick';
                SX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
                patchData.Vertices = obj(index).QN1;
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

        %% helpful methods
        %        function checkIntegerMatrix2(obj,input)
        %            assert(isnumeric(input) & (sum(sum(mod(input,1))) == 0),'input need to be of type integer even if it is classified as float');
        %        end
        %% solver initialization
        out = assignGlobalDofs(obj,input)
        initializeQV(obj,input)
        assignShapefunctions(obj)
        %% mass matrix routine
        out = massMatrix(obj)
        out = accessElements(obj,input)
    end
end