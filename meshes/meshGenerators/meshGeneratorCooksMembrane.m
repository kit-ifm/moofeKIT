function [NODES,EDOF,bounEDOFs] = meshGeneratorCooksMembrane(nelX, nelY, nelZ, order, serendipity)
%% Netzgenerator fuer Quader mit triilinearen isoparametrische 8-Knotenelementen
%[NODES,EDOF,bounEDOFs] = meshGeneratorCube(lengthX, lengthY, lengthZ, nelx, nely, nelz)
%
%INPUT
%--------------------------------------
% lengthX: length of brick in x-direction
% numberOfElementsX: number of Elements in x-direction
% order: interpolation order of the elements (may be 1 or 2)
% serendipity: if serendipity elements are to be used for order==2
%
%OUTPUT
%--------------------------------------
% NODES: positon of nodes
% EDOF: elements (trilinear)
% BOUNEDOF: structure with boundary elements (bilinear)
%   .SX1 surface perpendicular to x-axis at x=xmin
%   .SX2 surface perpendicular to x-axis at x=xmax


% F.STREICH
% Last changes: 27.04.2017 S.Schuss, 27.11.2017 R. Pfefferkorn
% 3.9.2020 R. Pfefferkorn: added support for quadratic elements
assert(floor(nelX)==nelX & nelX>0,'nelx must be an integer and >0')
assert(floor(nelY)==nelY & nelY>0,'nely must be an integer and >0')
assert(floor(nelZ)==nelZ & nelZ>0,'nelz must be an integer and >0')

%set standard arguments
if nargin<=3
    order = 1;
end
if nargin<=4
    serendipity = false;
end
if order==1 && serendipity==true
    warning('setting serendipity=true has no effect if order is 1')
end

%% NODES
%number of nodes
nel = nelX*nelY*nelZ;
nnoX = nelX*order+1;
nnoY = nelY*order+1;
nnoZ = nelZ*order+1;
%geometry
heightLeft = 44;
heightRight = 16;
length = 48;
thickness = 4;

geometry = [0 length; 0 thickness; 0 heightLeft];
nodesXDir = linspace(geometry(1,1),geometry(1,2),nnoX)';
nodesYDir = linspace(geometry(2,1),geometry(2,2),nnoY)';
nodesZDir = linspace(geometry(3,1),geometry(3,2),nnoZ)';
zw = [kron(ones(nnoY,1),nodesXDir), kron(nodesYDir,ones(nnoX,1))]; %nodes 2D
NODES = [kron(ones(nnoZ,1),zw), kron(nodesZDir,ones(size(zw,1),1))];

NODES(:,3) = NODES(:,3) + NODES(:,1)/length*heightLeft + (NODES(:,3)/heightLeft).*(NODES(:,1)/length*(heightRight-heightLeft));


%% EDOF
%e.dof 2D
if order == 1
    edof2D = meshEdofRectangle(1, nelX, nelY, order);

    %edof
    EDOF = zeros(nel,8);
    EDOF(1:nelX*nelY,:) = [edof2D, edof2D+nnoX*nnoY];
    EDOF(:,:) = kron(ones(nelZ,1),EDOF(1:nelX*nelY,:))+kron((0:nelZ-1)',order*nnoX*nnoY*ones(nelX*nelY,8));
elseif order == 2
    %2d mesh 9-node quadratic
    edof2D = meshEdofRectangle(1, nelX, nelY, order);
    
    %edof 3D unstructured for one layer (only 2d structuring is applied)
    edofZw = zeros(nelX*nelY,27);
    edofZw(1:nelX*nelY,:) = [edof2D, edof2D+nnoX*nnoY, edof2D+2*nnoX*nnoY];
    
    %structure first layer with node numbering of esra
    %1-8 corner nodes
    %9-16 midside nodes
    %17 midplane node (irregular)
    %18-21 midside nodes
    %22-26 midplane nodes
    %27 center node
    EDOF = zeros(nel,27);
    EDOF(1:nelX*nelY,:) = edofZw(:,[1,2,3,4,19,20,21,22, 10,11,12,13,5,6,7,8, 9, 23,24,25,26, 27,14,15,16,17, 18]);
    %complete EDOF
    EDOF(:,:) = kron(ones(nelZ,1),EDOF(1:nelX*nelY,:))+kron((0:nelZ-1)',order*nnoX*nnoY*ones(nelX*nelY,27));
else
    error('order not implemented')
end

%% boundary Edofs
% get the elements that touch the border, then select the nodes on the
% boundary in the correct order
if order == 1
    elementSet = kron(ones(nelZ,1),(1:nelX:nelX*nelY-nelX+1)')+kron((0:nelZ-1)',ones(nelY,1))*nelX*nelY;
    bounEDOFs.SX1 = EDOF(elementSet, [4,1,5,8]);
    bounEDOFs.SX2 = EDOF(elementSet+nelX-1, [2,3,7,6]);

    elementSet = kron(ones(nelZ,1),(1:nelX)')+kron((0:nelZ-1)',ones(nelX,1))*nelX*nelY;
    bounEDOFs.SY1 = EDOF(elementSet, [1,2,6,5]);
    bounEDOFs.SY2 = EDOF(elementSet+nelX*nelY-nelX, [3,4,8,7]);

    elementSet = 1:nelX*nelY;
    bounEDOFs.SZ1 = EDOF(1:nelX*nelY, [4,3,2,1]);
    bounEDOFs.SZ2 = EDOF(elementSet+nelX*nelY*nelZ-nelX*nelY, [5,6,7,8]);
elseif order == 2
    elementSet = kron(ones(nelZ,1),(1:nelX:nelX*nelY-nelX+1)')+kron((0:nelZ-1)',ones(nelY,1))*nelX*nelY;
    bounEDOFs.SX1 = EDOF(elementSet, [4,1,5,8, 16,9,21,12, 26]);
    bounEDOFs.SX2 = EDOF(elementSet+nelX-1, [2,3,7,6, 14,11,19,10, 24]);
    
    elementSet = kron(ones(nelZ,1),(1:nelX)')+kron((0:nelZ-1)',ones(nelX,1))*nelX*nelY;
    bounEDOFs.SY1 = EDOF(elementSet, [1,2,6,5, 13,10,18,9, 23]);
    bounEDOFs.SY2 = EDOF(elementSet+nelX*nelY-nelX, [3,4,8,7, 15,12,20,11, 25]);

    elementSet = 1:nelX*nelY;
    bounEDOFs.SZ1 = EDOF(1:nelX*nelY, [4,3,2,1, 15,14,13,16, 17]);
    bounEDOFs.SZ2 = EDOF(elementSet+nelX*nelY*nelZ-nelX*nelY, [5,6,7,8, 18,19,20,21, 22]);
else
    error('order not implemented')
end

%change mesh for serendipity
if order==2 && serendipity
    %reduce edofs
    EDOF = EDOF(:,[1:8, 13:16, 18:21, 9:12]);
    bounEDOFs.SX1 = bounEDOFs.SX1(:, 1:8);
    bounEDOFs.SX2 = bounEDOFs.SX2(:, 1:8);
    bounEDOFs.SY1 = bounEDOFs.SY1(:, 1:8);
    bounEDOFs.SY2 = bounEDOFs.SY2(:, 1:8);
    bounEDOFs.SZ1 = bounEDOFs.SZ1(:, 1:8);
    bounEDOFs.SZ2 = bounEDOFs.SZ2(:, 1:8);
    
    %find surplus nodes
    delNodes = (1:size(NODES,1))';
    delNodes(unique(EDOF)) = [];
    
    %delete surplus nodes
    [EDOF,NODES] = deleteNodes(EDOF,NODES,delNodes);
    bounEDOFs.SX1 = deleteNodes(bounEDOFs.SX1,[],delNodes);
    bounEDOFs.SX2 = deleteNodes(bounEDOFs.SX2,[],delNodes);
    bounEDOFs.SY1 = deleteNodes(bounEDOFs.SY1,[],delNodes);
    bounEDOFs.SY2 = deleteNodes(bounEDOFs.SY2,[],delNodes);
    bounEDOFs.SZ1 = deleteNodes(bounEDOFs.SZ1,[],delNodes);
    bounEDOFs.SZ2 = deleteNodes(bounEDOFs.SZ2,[],delNodes);
end


end

function [EDOF,NODES] = deleteNodes(EDOF,NODES,delNodes)
    %deletes all nodes from given in delNodes
    delNodes = unique(delNodes);

    %delete all the node
    if nargout==2
        NODES(delNodes,:) = [];
    end
    
    %shift the node count in the EDOF
    for i=1:size(delNodes,1)
        % delete node
        actDupl = delNodes(i);
        EDOF(EDOF>=actDupl) = EDOF(EDOF>=actDupl)-1;
        delNodes(i:end) = delNodes(i:end)-1;
    end
end