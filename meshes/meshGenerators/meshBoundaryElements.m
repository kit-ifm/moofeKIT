function [nodes,edof] = meshBoundaryElements(inpnodes, line, tol, order)
%gives node list and edof of boundary elements (works only if boundary is
%indexed ascendingly, edof wrong otherwise)
%   input: nodes (all nodes from which subset should be selected)
%   input: line (2D-array with points on which nodes are searched)
%   input: tol (tolerance for node search)
%   input: order (order of the FE-mesh)
%   output: nodes (all nodes on the boundary)
%   output: edof (edof of boundary elements)

%% nodes
nodes = findNodes(inpnodes,line,tol); 
nno = size(nodes,1);

%% edof
edof = zeros((nno-1)/order,order+1);
switch order
    case 1
        edof(:,1)=nodes(1:end-1);
        edof(:,2)=nodes(2:end);
    case 2
        edof(:,1)=nodes(1:2:end-2);
        edof(:,2)=nodes(3:2:end);
        edof(:,3)=nodes(2:2:end-1);
    otherwise
        error('order not impelemented')
end
end

