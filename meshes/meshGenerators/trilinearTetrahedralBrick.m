function [NODES,edofTet,edofneumannTet] = trilinearTetrahedralBrick(lengthX, lengthY, lengthZ, nelx, nely, nelz)

%%Startscript%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functionname:     trilinearTetrahedralBrick.m
% Creation date:    09.03.2018
% Creator:          P. Kinon 
% Discription:      3D-mesh-generation for 4-node trilinear tet-elements in a
%                   hexahedron geometry (e.g. cooks membrane).
%                   Edof, nodes and neumann-edof are created by using
%                   geometry data.
% Based on:         10-node trilinear hex-elements "brick3D.m"

% Modifications:    --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(lengthX>0,'lengthX must be positive')
assert(lengthY>0,'lengthY must be positive')
assert(lengthZ>0,'lengthZ must be positive')
assert(floor(nelx)==nelx & nelx>0,'nelx must be an integer and >0')
assert(floor(nely)==nely & nely>0,'nely must be an integer and >0')
assert(floor(nelz)==nelz & nelz>0,'nelz must be an integer and >0')

%% NODES
%number of nodes
nnoX = nelx+1;
nnoY = nely+1;
nnoZ = nelz+1;
%geometry
geometry = [-lengthX/2 lengthX/2; -lengthY/2 lengthY/2; -lengthZ/2 lengthZ/2];
nodesXDir = linspace(geometry(1,1),geometry(1,2),nnoX)';
nodesYDir = linspace(geometry(2,1),geometry(2,2),nnoY)';
nodesZDir = linspace(geometry(3,1),geometry(3,2),nnoZ)';
zw = [kron(ones(nnoY,1),nodesXDir), kron(nodesYDir,ones(nnoX,1))]; %nodes 2D
NODES = [kron(ones(nnoZ,1),zw), kron(nodesZDir,ones(size(zw,1),1))];

%% EDOF
%edof 2D
edof2D = zeros(nelx*nely,4);
edof2D(:,1)=kron(ones(nely,1),(1:(nnoX-1))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelx,1));
edof2D(:,2)=kron(ones(nely,1),(2:(nnoX-0))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelx,1));
edof2D(:,3)=kron(ones(nely,1),(2:(nnoX-0))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelx,1));
edof2D(:,4)=kron(ones(nely,1),(1:(nnoX-1))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelx,1));

%edof
EDOF = zeros(nelx*nely*nelz,8);
EDOF(1:nelx*nely,:) = [edof2D, edof2D+nnoX*nnoY];
EDOF(:,:) = kron(ones(nelz,1),EDOF(1:nelx*nely,:))+kron((0:nelz-1)',nnoX*nnoY*ones(nelx*nely,8));

%% boundary Edofs
elementSet = kron(ones(nelz,1),(1:nelx:nelx*nely-nelx+1)')+kron((0:nelz-1)',ones(nely,1))*nelx*nely;
bounEDOFs.SX1 = EDOF(elementSet, [4,1,5,8]);
bounEDOFs.SX2 = EDOF(elementSet+nelx-1, [2,3,7,6]);

elementSet = kron(ones(nelz,1),(1:nelx)')+kron((0:nelz-1)',ones(nelx,1))*nelx*nely;
bounEDOFs.SY1 = EDOF(elementSet, [1,2,6,5]);
bounEDOFs.SY2 = EDOF(elementSet+nelx*nely-nelx, [8,7,3,4]);

elementSet = 1:nelx*nely;
bounEDOFs.SZ1 = EDOF(1:nelx*nely, [4,3,2,1]);
bounEDOFs.SZ2 = EDOF(elementSet+nelx*nely*nelz-nelx*nely, [5,6,7,8]);

%--------------------------------------------------------------
% Einteilung eines Hexaeder-Elements in je 5 Tetraeder-Elemente
%--------------------------------------------------------------

newEl=size(EDOF,1)*5;
edofTet=zeros(newEl,4);
edofneumannTet=zeros(2*size(bounEDOFs.SZ2,1),3);

for p=1:(newEl/5)
    edofTet(5*(p-1)+1,:)=EDOF(p,[8 5 1 6]);
    edofTet(5*(p-1)+2,:)=EDOF(p,[8 7 6 3]);
    edofTet(5*(p-1)+3,:)=EDOF(p,[3 4 1 8]);
    edofTet(5*(p-1)+4,:)=EDOF(p,[3 2 6 1]);
    edofTet(5*(p-1)+5,:)=EDOF(p,[3 6 8 1]);
end

for pp=1:size(bounEDOFs.SZ2,1)
    edofneumannTet(2*(pp-1)+1,:)=bounEDOFs.SZ2(pp,[1 2 4]);
    edofneumannTet(2*(pp-1)+2,:)=bounEDOFs.SZ2(pp,[3 4 2]);
end

end