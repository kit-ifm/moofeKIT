function [NODES, EDOF, bounEDOFs] = meshRectangle(lengthX, lengthY, numberOfElementsX, numberOfElementsY, order, serendipity)
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
assert(lengthX > 0, 'lengthX must be positive')
assert(lengthY > 0, 'lengthY must be positive')
assert(floor(numberOfElementsX) == numberOfElementsX & numberOfElementsX > 0, 'numberOfElementsX must be an integer and >0')
assert(floor(numberOfElementsY) == numberOfElementsY & numberOfElementsY > 0, 'numberOfElementsY must be an integer and >0')

%set standard arguments
if nargin <= 4
    order = 1;
end
if nargin <= 5
    serendipity = false;
end
if order == 1 && serendipity == true
    warning('setting serendipity=true has no effect if order is 1')
end

%% Nodes
% number of nodes (total)
nelX = numberOfElementsX; % numberOfElements in x-direction
nelY = numberOfElementsY; % numberOfElements in y-direction
nnoX = numberOfElementsX * order + 1; % numberOfNodes in x-direction
nnoY = numberOfElementsY * order + 1; % numberOfNodes in y-direction

x = -lengthX / 2:lengthX / (numberOfElementsX * order):lengthX / 2;
y = -lengthY / 2:lengthY / (numberOfElementsY * order):lengthY / 2;
NODES = zeros(nnoX*nnoY, 2);

NODES(:, 1) = kron(ones(1, numberOfElementsY*order+1), x);
NODES(:, 2) = kron(y, ones(1, numberOfElementsX*order+1));

%% Edof
EDOF = meshEdofRectangle(1, numberOfElementsX, numberOfElementsY, order);

%% edofs of boundary
bounEDOFs = struct('SY1', zeros(nelX, order+1), 'SY2', zeros(nelX, order+1), ...
    'SX1', zeros(nelY, order+1), 'SX2', zeros(nelY, order+1));

for i=1:order + 1
    bounEDOFs.SY1(:, i) = i:order:nnoX - order - 1 + i;
    bounEDOFs.SY2(:, i) = nnoX * nnoY - nnoX + i:order:nnoX * nnoY - order - 1 + i;
    bounEDOFs.SX1(:, i) = nnoX * (i-1) + 1:order*nnoX:nnoX * nnoY - nnoX * (order+2-i) + 1;
    bounEDOFs.SX2(:, i) = nnoX + nnoX * (i-1):order * nnoX:nnoX * nnoY - nnoX * (order+1-i);
end

%change mesh for serendipity
if serendipity
    if order == 2
        % reduce edofs
        EDOF = EDOF(:, 1:8);

        % find surplus nodes
        delNodes = (1:size(NODES, 1))';
        delNodes(unique(EDOF)) = [];

        % delete surplus nodes
        [EDOF, NODES] = deleteNodes(EDOF, NODES, delNodes);
        bounEDOFs.SX1 = deleteNodes(bounEDOFs.SX1, [], delNodes);
        bounEDOFs.SX2 = deleteNodes(bounEDOFs.SX2, [], delNodes);
        bounEDOFs.SY1 = deleteNodes(bounEDOFs.SY1, [], delNodes);
        bounEDOFs.SY2 = deleteNodes(bounEDOFs.SY2, [], delNodes);
    elseif order > 2
        error(['Serendipity elements currently not implemented for order ', num2str(order), '!']);
    end
end

end

function [EDOF, NODES] = deleteNodes(EDOF, NODES, delNodes)
%deletes all nodes from given in delNodes
delNodes = unique(delNodes);

%delete all the node
if nargout == 2
    NODES(delNodes, :) = [];
end

%shift the node count in the EDOF
for i = 1:size(delNodes, 1)
    % delete node
    actDupl = delNodes(i);
    EDOF(EDOF >= actDupl) = EDOF(EDOF >= actDupl) - 1;
    delNodes(i:end) = delNodes(i:end) - 1;
end
end