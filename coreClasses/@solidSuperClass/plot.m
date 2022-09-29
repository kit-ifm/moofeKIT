function plot(obj,setupObject)
% PLOT function for ploting continuum classes.
%
% CALL
% plot(obj,setupObject,varargin)
% obj: The first argument is expected to be a plotable object like
%       instanciated from class continuumClass.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
%
% CREATOR(S)
% Christian Hesch, Jonathan Schulte (esra) 11.2015; Marlon Franke
% FIXME: needs to be rebuild

for i = 1:numel(obj)
    switch setupObject.plotObject.time
        case 'R'
            q = obj.qR;
        case 'N'
            q = obj.qN;
        case 'N1'
            q = obj.qN1;
    end
    if isa(obj, 'thermoClass')
        ed = obj.meshObject.nodes;
    else
        ed = q;
    end
    obj(i).plotData.detColorDataFcn = @detColorData;
    colorData = obj(i).plotData.detColorDataFcn(obj(i),setupObject, q);
    [patchData, meshData, patchFlag] = createPatchLagrange(obj(i),setupObject,colorData, ed);
    [patchData, meshData] = slimPatch(patchData, meshData);        % remove double vertices and patches
    if patchFlag
        patch(patchData);
    end
    colorbar;
    %% scale caxis
    if ~strcmpi(setupObject.plotObject.colorData,'none') && isempty(setupObject.plotObject.colorBarLimits)
        clim = caxis;
        clim = [min(min(colorData,clim(1))),max(max(colorData,clim(2)))];
        clim = [min(clim(1),0),max(clim(2),0)];
        if clim(1)~=clim(2)
            caxis(clim)
        end
    end    
end
end
%% Function to determine color data
function colorData = detColorData(obj,setupObject, q)
colorData = [];
plotObject = setupObject.plotObject;
if ~strcmpi(plotObject.postPlotType,'none')
    fdof = obj.meshObject.globalNodesDof(:,1);
    maxElAct = max(fdof);
    if strcmpi(plotObject.postPlotType,'stress') || strcmpi(plotObject.postPlotType,'D1') || strcmpi(plotObject.postPlotType,'D2') || strcmpi(plotObject.postPlotType,'D3')
        postDataFE = callElements(obj, setupObject, true);        
        R = sparse(vertcat(postDataFE(:).indexReI),1,vertcat(postDataFE(:).Se),maxElAct,1);
        K = sparse(vertcat(postDataFE(:).indexKeI), vertcat(postDataFE(:).indexKeJ), vertcat(postDataFE(:).Me),maxElAct,maxElAct);
        colorData = zeros(size(fdof,1),1);
        colorData(:) = K(fdof,fdof)\R(fdof,1);
    else
        factor = 1;
        fieldPosition = 1;
        if ~isa(obj,'thermoClass')
            fieldPosition = obj.dimension;
        end
        if strcmpi(plotObject.postPlotType,'phi')
            fieldPosition = fieldPosition + 1;
        elseif strcmpi(plotObject.postPlotType,'temp')
            if isa(obj,'solidThermoClass')
                fieldPosition = fieldPosition + 1;
            elseif isa(obj,'solidElectroThermoClass')
                fieldPosition = fieldPosition + 2;
            end
        elseif strcmpi(plotObject.postPlotType,'zero')
            factor = 0;
        end
        colorData = q(:,fieldPosition)*factor;
    end
end
end