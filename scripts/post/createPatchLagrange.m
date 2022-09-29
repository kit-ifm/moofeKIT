function [patchData, meshData, patchFlag] = createPatchLagrange(obj,setupObject,colorData, ed)
edof = obj.meshObject.edof;
dimension = obj.dimension;
patchData = struct('Vertices',[],'Faces',[]);
patchFlag = true;
%% detect element type
if size(edof,2) == 8 && dimension ==3 || size(edof,2) == 27 && dimension == 3
    element = '3DBrick';
    SX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
elseif size(edof,2) == 20 && dimension == 3
    element = '3DSerendipity';
    SX = [1 9 2 18 6 13 5 17;
        2 10 3 19 7 14 6 18;
        3 11 4 20 8 15 7 19;
        4 12 1 17 5 16 8 20 ;
        5 13 6 14 7 15 8 16 ;
        1 9 2 10 3 11 4 12 ];
% elseif size(edof,2) == 4 && dimension == 3 && isSurface %Surface of 3D Brick
%     element = '2DQuadratic';
%     SX = [1 2 3 4];
elseif size(edof,2) == 4 && dimension == 3
    element = '3DTet';
    SX = [1 2 3; 1 2 4; 2 3 4; 1 3 4;];
elseif size(edof,2) == 10 && dimension == 3
    element = '3DTetQuadratic';
    SX = [1 5 2 6 3 7; 1 5 2 9 4 8; 2 6 3 10 4 9; 1 7 3 10 4 8];
elseif size(edof,2) == 6
    element = '2DTetQuadratic';
     SX = [1 6 3 5 2 4];
elseif size(edof,2) == 4 && dimension == 2
    element = '2DBrick';
    SX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
elseif size(edof,2) == 8 && dimension == 2
    element = '2DSerendipity';
    SX = [1 5 2 6 3 7 4 8];
elseif size(edof,2) == 9 && dimension ==2 %|| (dimension == 3 && isSurface)) % 2D - Brick
    element = '2DQuadratic';
    SX = [1 2 3 4 5 6 7 8 9];
elseif size(edof,2) == 3 && dimension == 2
    element = '2DTet';
    SX = [1 2 3; 1 2 4; 2 3 4; 1 3 4;];
elseif dimension == 1
    element = '1D';
else
    error('Element-type not implemented')
end
%%Faces
switch element
    case {'3DBrick','3DSerendipity'}
        patchData.Faces = [edof(:,SX(1,:));
            edof(:,SX(2,:));
            edof(:,SX(3,:));
            edof(:,SX(4,:));
            edof(:,SX(5,:));
            edof(:,SX(6,:))];
    case {'3DTet','3DTetQuadratic'}
        patchData.Faces = [edof(:,SX(1,:));
            edof(:,SX(2,:));
            edof(:,SX(3,:));
            edof(:,SX(4,:));];
    case {'2DBrick', '2DTet','2DSerendipity','2DTetQuadratic'}
        patchData.Faces = edof(:,SX(1,:));
    case {'2DQuadratic'}
        patchData.Faces = edof(:,SX(1,1:4));
    case '1D'
        numberOfElements = size(edof, 1);
        numberOfElementNodes = size(edof, 2);
        edgeColor = 'interp';
        faceStyle = 'none';
        % plot line
        x = zeros(numberOfElements, numberOfElementNodes);
        c = zeros(numberOfElements, numberOfElementNodes);
        x(:) = ed(edof);
        y = zeros(size(x));
        c(:) = colorData(edof);
        patch(x', y', c', 'edgecolor', edgeColor, 'linewidth', 2*setupObject.plotObject.lineWeight, 'facecolor', faceStyle, 'facealpha', 0.9);
        % plot outside nodes
        xOutsideNodes = zeros(numberOfElements, 2);
        cOutsideNodes = zeros(numberOfElements, 2);
        xOutsideNodes(:) = ed(edof(:, 1:2));
        yOutsideNodes = zeros(size(xOutsideNodes));
        cOutsideNodes(:) = colorData(edof(:, 1:2));
        patch(xOutsideNodes', yOutsideNodes', cOutsideNodes', 'edgecolor', edgeColor, 'marker', '.', 'markersize', 15, 'facecolor', faceStyle, 'facealpha', 0.9, 'LineStyle', 'none');
        % plot inner nodes
        if numberOfElementNodes > 2
            xInnerNodes = zeros(numberOfElements, numberOfElementNodes-2);
            cInnerNodes = zeros(numberOfElements, numberOfElementNodes-2);
            xInnerNodes(:) = ed(edof(:, 3:end));
            yInnerNodes = zeros(size(xInnerNodes));
            cInnerNodes(:) = colorData(edof(:, 3:end));
            patch(xInnerNodes', yInnerNodes', cInnerNodes', 'edgecolor', edgeColor, 'marker', '.', 'markersize', 8, 'facecolor', faceStyle, 'facealpha', 0.9, 'LineStyle', 'none');
        end
        patchFlag = false;
end
switch setupObject.plotObject.time
    case {'R', 'N', 'N1'}
        patchData.Vertices = ed(:,1:dimension);
    otherwise
        assert(isscalar(options.time),'Please give timestep with scalar')
        timestep = options.time;
end
patchData.FaceVertexCData = colorData;
%     patchData.FaceColor = [0 1 0];
patchData.FaceAlpha = 0.5;
patchData.FaceColor = 'interp';
patchData.EdgeColor = 'black';
meshData = [];
end