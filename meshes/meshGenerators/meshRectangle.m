function [NODES,EDOF,bounEDOFs] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order)
% mesh for a rectangle with dimensions [-lengthX/2,lengthX/2] x [-lengthY/2,lengthY/2]

%INPUT
%--------------------------------------
% lengthX: length of brick in x-direction
% numberOfElementsX: number of Elements in x-direction
% order: interpolation order of the elements (may be 1 or 2)

%OUTPUT
%--------------------------------------
% NODES: positon of nodes
% EDOF: elements (trilinear)
% BOUNEDOF: structure with boundary elements (bilinear)
%   .SX1 surface perpendicular to x-axis at x=xmin
%   .SX2 surface perpendicular to x-axis at x=xmax

% Robin Pfefferkorn 13.6.17
% Last Change: 26.11.17 Robin Pfefferkorn (added edofs of boundarys)


%% check input
assert(lengthX>0,'lengthX must be positive')
assert(lengthY>0,'lengthY must be positive')
assert(floor(numberOfElementsX)==numberOfElementsX & numberOfElementsX>0,'numberOfElementsX must be an integer and >0')
assert(floor(numberOfElementsY)==numberOfElementsY & numberOfElementsY>0,'numberOfElementsY must be an integer and >0')
assert(order==1 | order==2,'order of mesh must be 1 or 2 (bilinear or biquadratic)')

%% Nodes
% number of nodes (total)
nelX = numberOfElementsX;
nelY = numberOfElementsY;
nnoX=numberOfElementsX*order+1; %numberOfNodes
nnoY=numberOfElementsY*order+1; %numberOfNodes

x = -lengthX/2:lengthX/(numberOfElementsX*order):lengthX/2;
y = -lengthY/2:lengthY/(numberOfElementsY*order):lengthY/2;
NODES = zeros(nnoX*nnoY,2);

NODES(:,1) = kron(ones(1,numberOfElementsY*order+1),x);
NODES(:,2) = kron(y,ones(1,numberOfElementsX*order+1));

%% Edof
EDOF = meshEdofRectangle(1,numberOfElementsX,numberOfElementsY,order);

%% edofs of boundary
bounEDOFs = struct('SY1',zeros(nelX,order+1),'SY2',zeros(nelX,order+1),...
    'SX1',zeros(nelY,order+1),'SX2',zeros(nelY,order+1));

if order == 1 %bilinear shape functions
    bounEDOFs.SY1(:,1) = 1:1:nnoX-1;
    bounEDOFs.SY1(:,2) = 2:1:nnoX-0;
    
    bounEDOFs.SY2(:,1) = nnoX*nnoY-nnoX+1:1:nnoX*nnoY-1;
    bounEDOFs.SY2(:,2) = nnoX*nnoY-nnoX+2:1:nnoX*nnoY-0;
    
    bounEDOFs.SX1(:,1) = nnoX*0+1:nnoX:nnoX*nnoY-nnoX*2+1;
    bounEDOFs.SX1(:,2) = nnoX*1+1:nnoX:nnoX*nnoY-nnoX*1+1;
    
    bounEDOFs.SX2(:,1) = nnoX+nnoX*0:nnoX:nnoX*nnoY-nnoX*1;
    bounEDOFs.SX2(:,2) = nnoX+nnoX*1:nnoX:nnoX*nnoY-nnoX*0;
elseif order == 2 %biquadratic shape functions
    bounEDOFs.SY1(:,1) = 1 : 2 : nnoX-2;
    bounEDOFs.SY1(:,2) = 2 : 2 : nnoX-1;
    bounEDOFs.SY1(:,3) = 3 : 2 : nnoX-0;
    
    bounEDOFs.SY2(:,1) = nnoX*nnoY-nnoX+1 : 2 : nnoX*nnoY-2;
    bounEDOFs.SY2(:,2) = nnoX*nnoY-nnoX+2 : 2 : nnoX*nnoY-1;
    bounEDOFs.SY2(:,3) = nnoX*nnoY-nnoX+3 : 2 : nnoX*nnoY-0;
    
    bounEDOFs.SX1(:,1) = nnoX*0+1 : 2*nnoX : nnoX*nnoY-nnoX*3+1;
    bounEDOFs.SX1(:,2) = nnoX*1+1 : 2*nnoX : nnoX*nnoY-nnoX*2+1;
    bounEDOFs.SX1(:,3) = nnoX*2+1 : 2*nnoX : nnoX*nnoY-nnoX*1+1;
    
    bounEDOFs.SX2(:,1) = nnoX+nnoX*0 : 2*nnoX : nnoX*nnoY-nnoX*2;
    bounEDOFs.SX2(:,2) = nnoX+nnoX*1 : 2*nnoX : nnoX*nnoY-nnoX*1;
    bounEDOFs.SX2(:,3) = nnoX+nnoX*2 : 2*nnoX : nnoX*nnoY-nnoX*0;
else
    error('unable to mesh: elementtype (order) currently not implemented') 
end

end