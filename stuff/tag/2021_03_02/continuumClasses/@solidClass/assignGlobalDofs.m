function out = assignGlobalDofs(obj,input)
%% global dofs
% general values
for index = 1:numel(obj)
    dimension = obj(index).dimension;
    edof = obj(index).edof;
    numberOfNodes = size(obj(index).nodes,1);
    numberOfElements = size(edof,1);
    globalFullEdof = zeros(numberOfElements,size(edof,2)*dimension);
    for j = 1:dimension
        globalFullEdof(:,j:dimension:end) = input + edof*dimension-(dimension-j);
    end
    globalNodesDOF = zeros(numberOfNodes,dimension);
    globalNodesDOF(:,1:dimension) = input + kron(ones(numberOfNodes,1),1:dimension)+kron((0:dimension:dimension*numberOfNodes-1)',ones(1,dimension));
    obj(index).globalFullEdof = globalFullEdof;
    obj(index).globalNodesDOF = globalNodesDOF;    
    input = input + numberOfNodes*dimension;
end
out = input;
end
