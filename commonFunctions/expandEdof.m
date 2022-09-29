function [edof1,edof2] = expandEdof(edof)
edofDouble = double(edof);
edof1 = kron(ones(1,size(edofDouble,2)),edofDouble);
edof2 = kron(edofDouble,ones(1,size(edofDouble,2)));
end