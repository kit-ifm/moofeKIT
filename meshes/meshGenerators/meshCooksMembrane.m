function [nodes, edof, edofNeumann] = meshCooksMembrane(numberOfElementsX, numberOfElementsY, order, serendipity)
%MESHCOOKSMEMBRANE Mesh for the Cooks Membrane
%   This function returns the nodes and the edof for the Cooks Membrane
%
%   CALL
%   [nodes, edof, bounEdof] = meshCooksMembrane(numberOfElementsX, numberOfElementsY, order, serendipity)
%   numberOfElementsX: number of elements in X-direction
%   numberOfElementsY: number of elements in Y-direction
%   order: order of the ansatz functions
%   serendipity: serendipity shape functions (true / false)
%   nodes: coordinates of the nodes
%   edof: node numbers for all elements
%   bounEdof: Neumann boundary (right side of membrane)
%
%   REFERENCE
%   -
%
%   CREATOR(S)
%   Jakob Hammes

% dimensions
heightLeft = 44;
heightRight = 16;
length = 48;

% create rectangle mesh first 
[nodesOriginal, edof, bounEdof] = meshRectangle(length, heightLeft, numberOfElementsX, numberOfElementsY, order, serendipity);
nodesOriginal(:,1) = nodesOriginal(:,1) + length/2; 
nodesOriginal(:,2) = nodesOriginal(:,2) + heightLeft/2;
edofNeumann = bounEdof.SX2;

% distort rectangle mesh
nodes = nodesOriginal;
nodes(:,2) = nodesOriginal(:,2) + nodesOriginal(:,1)/length*heightLeft + (nodesOriginal(:,2)/heightLeft).*(nodesOriginal(:,1)/length*(heightRight-heightLeft));
end