classdef solidThermoClass < solidSuperClass  
    properties (Constant = true)
        flagNewton = struct('initializeGlobalDofs',true,...
                            'initializeQV',true,...
                            'initializeShapeFunctions',true,...
                            'updateTimeDependentField',false,...
                            'updateContinuumFieldPreNewtonLoop',true,...
                            'callElements',true,...
                            'updateContinuumFieldNewtonLoop',true);
        additionalFields = 1;
    end
    properties
        elementDisplacementType = 'displacement';
    end
    methods
        %% constructor
        function obj = solidThermoClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
        function plot(obj)
            for index = 1:numel(obj)
                patchData = struct('Vertices',[],'Faces',[]);
                %     element = '3Dbrick';
                SX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
                patchData.Vertices = obj(index).qN1(:,1:obj(index).dimension);
                patchData.Faces = [ obj(index).edof(:,SX(1,:));
                    obj(index).edof(:,SX(2,:));
                    obj(index).edof(:,SX(3,:));
                    obj(index).edof(:,SX(4,:));
                    obj(index).edof(:,SX(5,:));
                    obj(index).edof(:,SX(6,:))];
                patchData.FaceVertexCData = obj(index).qN1(:,obj(index).dimension+1);
                patchData.FaceColor = 'interp';
                patchData.Linewidth = 0.5;
                patchData.EdgeColor = 'black';
                patchData.LineStyle = '-';
                patchData.FaceAlpha = 0.5;
                patch(patchData);
                hold on
            end
        end
    end
end