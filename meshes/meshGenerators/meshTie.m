function [nodes,edof] = meshTie(nodes1,edof1,nodes2,edof2,tolerance)
%MESHTIE ties mesh parts (deletes duplicate nodes)
%   tol...tolerance for deleting

%% check input
if size(nodes1,2) ~= size(nodes2,2)
    error('nodes1 and nodes2 must have the same dimensionality')
end
if size(edof1,2) ~= size(edof2,2)
    error('edof1 and edof2 must contain the same element types (same number of nodes)')
end

nodes = vertcat(nodes1,nodes2);
edof = vertcat(edof1,edof2);

nno = size(nodes,1);
dimension = size(nodes,2);

edofH = edof(:);
ind = 1:nno;

%% Find duplicates
x = zeros(nno,dimension);
distance = zeros(nno,1);

%array that contains replacement nodes, if zero no replacement
replaceNode = zeros(nno,1);

for i=1:size(nodes,1)-1
    %current node and distances to that node
    x(i:nno,:) = kron(ones(nno-i+1,1),nodes(i,:));
    distance(i) = 100*tolerance; %same node and nodes before can't be duplicates
    distance(i+1:nno) = sum(sqrt((nodes(i+1:nno,:)-x(i+1:nno,:)).^2),2);
    
    %find duplicates
    duplInd = ind(distance<=tolerance & replaceNode==0);
    
    if ~isempty(duplInd)
        replaceNode(duplInd) = i;
    end
end

duplicates = ind(replaceNode~=0);
replaceNode = replaceNode(duplicates);

%% Delete duplicate Nodes
% nodes(duplicates,:) = [];

for i=1:size(duplicates,2)
    % delete node
    actDupl = duplicates(i);
    nodes(actDupl,:) = [];
    
    % change node in edof
    edof(edof==actDupl) = replaceNode(i);
    
    % nodes in edof higher than actDupl have to be changed
    edof(edof>actDupl) = edof(edof>actDupl)-1;
    duplicates(duplicates>actDupl) = duplicates(duplicates>actDupl)-1;
    replaceNode(replaceNode>actDupl) = replaceNode(replaceNode>actDupl)-1;
end




