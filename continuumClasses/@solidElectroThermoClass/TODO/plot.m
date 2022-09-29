function plot(obj,varargin)
% Plot Object

% set function-handle for determination of the color data
for i = 1:numel(obj)
    obj(i).PLOT_DATA.detColorDataFcn = @detColorData;
end

plotObject(obj,varargin) %common/post/plotObject

end

%% Function to determine color data
function colorData = detColorData(obj,options)
%obj: your Object
%options: structure with all plot-options (see common/plot/ploObject/checkInputArguments.m)
colorData = [];
if options.phi
    switch options.time
        case 'REF'
            colorData =  obj.QREF(:,obj.DIM+1);
        case 'N'
            colorData =  obj.QN(:,obj.DIM+1);
        case 'N1'
            colorData =  obj.QN1(:,obj.DIM+1);
    end
elseif options.temp
    switch options.time
        case 'REF'
            colorData =  obj.QREF(:,obj.DIM+2);
        case 'N'
            colorData =  obj.QN(:,obj.DIM+2);
        case 'N1'
            colorData =  obj.QN1(:,obj.DIM+2);
    end
elseif options.stress~=0 || options.d1 || options.d2 || options.d3
    if strcmp(obj.DESIGNATOR,'mmhwSC')
        fdof = obj.NODESDOF(:,1);
        maxElAct = max(fdof);
        if options.stress~=0
            if options.d1
                internal=residual(obj,'d1',options.d1,'time',0);
            elseif options.d2
                internal=residual(obj,'d2',options.d2,'time',0);
            elseif options.d3
                internal=residual(obj,'d3',options.d3,'time',0);
            else
                internal=residual(obj,'stress',options.stress,'time',0);                
            end
        end
        rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
        kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
        colorData = zeros(size(fdof,1),1);
        colorData(:) = kHam(fdof,fdof)\rHam(fdof,1);
    end
end
end