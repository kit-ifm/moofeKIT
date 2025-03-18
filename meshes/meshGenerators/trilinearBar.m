function [nodes,edof,edofneumann] = trilinearBar(numberElementsX, numberElementsY, numberElementsZ, lenght, height, thickness)
%% net generator for trilinear, isoparametric 8-node elements (H1)
%   edof - element times dof-number in element e
%                   _ (dof number in element)  _
%           edof = | 1  4  5  2  10  13  14  11 | 1
%                  |            ...             | 2 (element)
%                  |            ...             | 3
%                  |_           ...            _| ...
%
%   nodes - node number times dimension
%                   _(x   y   z)_
%          nodes = |  0   0   0  | 1
%                  |  0  0.5  0  | 2 (node number)
%                  |     ...     | 3   
%                  |_    ...    _| ...
%           
% CREATOR(S)
% Bea Hummel
% 
% UPDATE
% Moritz Hille
%

%% generate nodes
dimension = 3;
numberNodes2D = (numberElementsX+1)*(numberElementsY+1); % number of nodes in profile (2D)
nodes = zeros(numberNodes2D*(numberElementsZ+1),dimension);
dx = lenght/numberElementsX;
dy = height/numberElementsY;
dz = thickness/numberElementsZ;

for i=1:numberElementsX+1
   x=(i-1)*dx;
   
   nodes((i-1)*(numberElementsY+1)+1:i*(numberElementsY+1),1)=x;                        %x
   nodes((i-1)*(numberElementsY+1)+1:i*(numberElementsY+1),2)=(0:numberElementsY)*dy;   %y
end

for i=1:numberElementsZ
    %x-y-Koordinate kopieren
    nodes(numberNodes2D*i+1:numberNodes2D*(i+1),1:2)=nodes(1:numberNodes2D,1:2);
    %z-Koordinate
    nodes(numberNodes2D*i+1:numberNodes2D*(i+1),3)=i*dz;
end

%% generate edof
edof = zeros(numberElementsX*numberElementsY*numberElementsZ,2^dimension);

for i=1:numberElementsX
   edof((i-1)*numberElementsY+1:i*numberElementsY,1:4)= ...
       [(1:numberElementsY)'+(i-1)*(numberElementsY+1),(1:numberElementsY)'+i*(numberElementsY+1), ...
        (1:numberElementsY)'+i*(numberElementsY+1)+1,(1:numberElementsY)'+(i-1)*(numberElementsY+1)+1];
end

edof(1:numberElementsX*numberElementsY,5:8) = edof(1:numberElementsX*numberElementsY,1:4)+numberNodes2D;
for i=2:numberElementsZ
    edof((i-1)*numberElementsX*numberElementsY+1:i*numberElementsX*numberElementsY,:) = edof(1:numberElementsX*numberElementsY,:)+numberNodes2D*(i-1);
end

%% generate edofneumann
edofneumann=zeros(numberElementsY*numberElementsZ,(dimension-1)^2);

for i=1:numberElementsZ
    edofneumann((i-1)*numberElementsY+1:i*numberElementsY,:)= ...
        [(1:numberElementsY)'+(numberElementsY+1)*numberElementsX+(i-1)*numberNodes2D, ...
         (1:numberElementsY)'+(numberElementsY+1)*numberElementsX+(i-1)*numberNodes2D+1, ...
         (1:numberElementsY)'+(numberElementsY+1)*numberElementsX+(i)*numberNodes2D+1, ...
         (1:numberElementsY)'+(numberElementsY+1)*numberElementsX+(i)*numberNodes2D];  
end