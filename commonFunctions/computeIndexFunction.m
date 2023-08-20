function [nodesXi, nodeIndeces] = computeIndexFunction(numberOfNodes) % checked for values numberOfNodes=4,9,16,25,36,49,64
% works only for SQUARE numberOfNodes, 2D elements
nodesXi = zeros(2, numberOfNodes);
nodeIndeces = zeros(numberOfNodes, 4); % array saves I-Y position per node + corresponding xi,eta coordinates
n = sqrt(numberOfNodes); % "size lenght" of element
UVec = uVec(n); % number of boundary nodes
% nxn elements
dim = n; % dimension of the nodes -> if node > U -> dim->dim-2;
counter = 0; % needed in order to switch dimensions

% The idea is the following: In most cases, not every node is situated
% on the boundary of the element. The nodes on the boundary are
% counted. This number is saved and the dimension is reduced.
% Virutally, the nodes on the boundary of the previous subdimension are now deleted and in a new, 'local'
% I,J coordinate system, starting with I=J=1, the nodes on the the new
% boundary are counted. And so on. This means a transformation has to
% take place for each 'subsdimension'.
for node = 1:numberOfNodes % loop over all nodes - > counting anti-clockwise on the boundary: question: How many nodes are there for each subdimension?
    [uLocal, nodeLocal, dim, counter] = nodesLocal(UVec, node, dim, counter); % transformation of the nodes in a local coordinate system (I,J)
    [dim] = nextDimension(node, dim, UVec, n);
    [I, J] = position(nodeLocal, dim); % position of the nodes in the local coordinate system
    [I, J] = localToGlobalCoordinates(I, J, n, dim); % transformation of the nodes in the global (I,J) coordinate system
    %disp([node,I,J])
    nodeIndeces(node, 1) = I; % saving coordinates
    nodeIndeces(node, 2) = J;

    xi = refCoordinate(I, n); % Calculating the corresponding xi,eta values
    eta = refCoordinate(J, n);
    nodesXi(1, node) = xi;
    nodesXi(2, node) = eta;
    nodeIndeces(node, 3) = xi;
    nodeIndeces(node, 4) = eta;
end


end


function [q, r] = euclideanDivision(n, d)

r = mod(n, d);
q = (n - r) / d;
end

function [I, J] = position(node, dim) % see 'Modular arithmetic'
[a, r] = euclideanDivision(node, 4); % gets the difference in position for every 'turn'; a a 'turn' means counting 4 nodes counter-clockwise

if mod(node, 4) == 1 % trivial properties
    I = 1 + a;
    J = 1;
elseif mod(node, 4) == 2
    I = dim;
    J = 1 + a;
elseif mod(node, 4) == 3
    I = dim - a;
    J = dim;
elseif mod(node, 4) == 0
    I = 1;
    J = dim + 1 - a;
end
end


function [counter] = numberOfU(n) %  computs the dimension of uVec, meaning the number of subdimensions
% (the Algorithm only considers nodes on the boundary, then virtually
% reduces the dimension meaning inner nodes are now on the boundary -> this
% step is repeated until all nodes have been considered.
dim = n;
counter = 0;
U = dim + 2 * (dim - 1) + (dim - 2);
while U > 0
    dim = dim - 2;
    U = dim + 2 * (dim - 1) + (dim - 2);
    counter = counter + 1;
end
if dim == 1
    counter = counter + 1;
end
end


function [uVec] = uVec(n) % boundary nodes for every subdimension
dim = numberOfU(n);
uVec = zeros(1, dim);
for ii = 1:dim
    U = n + 2 * (n - 1) + (n - 2); % in analogy to the circumference (german: "Umfang") for every subdimension
    % uVec contains the number of nodes on the boundary for every subdimension
    if U == 0 && mod(n, 2) == 1
        U = 1;
    end
    uVec(1, ii) = U;
    n = n - 2;
end
end

function [I, J] = localToGlobalCoordinates(I, J, n, dim) % local coordinates means: starting with I=J=1 in a subdimension
shift = ceil(0.5*(n - dim));
I = I + shift;
J = J + shift;
end


function [uLocal, nodeLocal, dim, counter] = nodesLocal(uVec, node, dim, counter) % calculates the local node number
Usum = uVec(1, 1);
if node <= uVec(1, 1)
    uLocal = uVec(1, 1);
    nodeLocal = node;
else
    uVecSize = size(uVec);
    for ii = 2:uVecSize(1, 2)
        Usum_old = Usum;
        Usum = Usum + uVec(1, ii);
        if node > Usum_old && node <= Usum
            if node == Usum %&& ii==uVecSize(1, 2)
                counter = 0;
            end
            nodeLocal = node - Usum_old;
            uLocal = uVec(ii);
            if counter == 0
                counter = counter + 1;
            end


        end
    end
end

end


function xi = refCoordinate(nodeIndex, numberOfNodes1D)
xi = (2 * nodeIndex - numberOfNodes1D - 1) / (numberOfNodes1D - 1);
end


function [dim] = nextDimension(node, dim, UVec, n)
uVecSize = size(UVec);
sum = 0;
for ii = 1:uVecSize(1, 2)
    sum = sum + UVec(ii);
    if node <= sum
        dim = n - (2 * (ii - 1));
        break
    end
end

end
