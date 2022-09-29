function [edof] = meshEdofRectangle(node1, nelX, nelY, order)

%% number elements and nodes (total)
nel = nelX * nelY; % total number of elements
nnoX = nelX * order + 1; % total number of nodes in x-direction

%% edof
if order == 1 % bilinear shape functions
    edof = zeros(nel, 4);
    edof(:, 1) = kron(ones(nelY, 1), (1:(nnoX - 1))') + kron((0:nnoX:(nnoX) * (nelY - 1))' + 0*nnoX, ones(nelX, 1));
    edof(:, 2) = kron(ones(nelY, 1), (2:(nnoX - 0))') + kron((0:nnoX:(nnoX) * (nelY - 1))' + 0*nnoX, ones(nelX, 1));
    edof(:, 3) = kron(ones(nelY, 1), (2:(nnoX - 0))') + kron((0:nnoX:(nnoX) * (nelY - 1))' + 1*nnoX, ones(nelX, 1));
    edof(:, 4) = kron(ones(nelY, 1), (1:(nnoX - 1))') + kron((0:nnoX:(nnoX) * (nelY - 1))' + 1*nnoX, ones(nelX, 1));
elseif order == 2 % biquadratic shape functions
    edof = zeros(nel, 9);
    edof(:, 1) = kron(ones(nelY, 1), (1:2:(nnoX - 2))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 0*nnoX, ones(nelX, 1));
    edof(:, 2) = kron(ones(nelY, 1), (3:2:(nnoX - 0))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 0*nnoX, ones(nelX, 1));
    edof(:, 3) = kron(ones(nelY, 1), (3:2:(nnoX - 0))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 2*nnoX, ones(nelX, 1));
    edof(:, 4) = kron(ones(nelY, 1), (1:2:(nnoX - 2))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 2*nnoX, ones(nelX, 1));
    edof(:, 5) = kron(ones(nelY, 1), (2:2:(nnoX - 1))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 0*nnoX, ones(nelX, 1));
    edof(:, 6) = kron(ones(nelY, 1), (3:2:(nnoX - 0))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 1*nnoX, ones(nelX, 1));
    edof(:, 7) = kron(ones(nelY, 1), (2:2:(nnoX - 1))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 2*nnoX, ones(nelX, 1));
    edof(:, 8) = kron(ones(nelY, 1), (1:2:(nnoX - 2))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 1*nnoX, ones(nelX, 1));
    edof(:, 9) = kron(ones(nelY, 1), (2:2:(nnoX - 1))') + kron((0:nnoX * 2:nnoX * 2 * (nelY - 1))' + 1*nnoX, ones(nelX, 1));
else
    error('unable to mesh: elementtype currently not implemented')
end

edof = edof + node1 - 1;
end