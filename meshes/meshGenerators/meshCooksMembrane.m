function [nodes,edof,edofNeumann] = meshCooksMembrane(numberOfElementsX, numberOfElementsY, order)
% Robin Pfefferkorn 30.4.19
% mesh for cooks membrane with dimensions:
heightLeft = 44;
heightRight = 16;
length = 48;

% create rectangle first 
x = 0:length/(numberOfElementsX*order):length;
y = 0:heightLeft/(numberOfElementsY*order):heightLeft;
nodes = zeros((numberOfElementsX*order+1)*(numberOfElementsY*order+1),2);

nodes(:,1) = kron(ones(1,numberOfElementsY*order+1),x);
nodes(:,2) = kron(y,ones(1,numberOfElementsX*order+1));
edof = meshEdofRectangle(1,numberOfElementsX,numberOfElementsY,order);

% distort rectangle mesh
nodes(:,2) = nodes(:,2) + nodes(:,1)/length*heightLeft + (nodes(:,2)/heightLeft).*(nodes(:,1)/length*(heightRight-heightLeft));

% edof Neumann
[~,edofNeumann] = meshBoundaryElements(nodes,[length,heightLeft;length,heightLeft+heightRight],length/numberOfElementsX*1e-5,order);
end