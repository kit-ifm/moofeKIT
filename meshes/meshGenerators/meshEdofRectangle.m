function edof = meshEdofRectangle(node1, nelX, nelY, order)

%% number elements and nodes (total)
nel = nelX * nelY; % total number of elements
nnoX = nelX * order + 1; % total number of nodes in x-direction

numberOfNodesPerElement = (order + 1)^2;

%% edof
edof = zeros(nel, numberOfNodesPerElement);
nodesXi = computeIndexFunction(numberOfNodesPerElement);

for ii = 1:numberOfNodesPerElement
    posX = round(order / 2 * (nodesXi(1, ii) + 1) + 1);
    posY = round(order / 2 * (nodesXi(2, ii) + 1) + 1);
    posXInv = round(order + 1 - posX);
    edof(:, ii) = kron(ones(nelY, 1), (posX:order:(nnoX - posXInv))') + kron((0:nnoX * order:nnoX * order * (nelY - 1))'+(posY - 1)*nnoX, ones(nelX, 1));
end

edof = edof + node1 - 1;
end