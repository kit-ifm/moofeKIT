function plot(obj,setupObject)
%plots boundaries from dirichlet

plotObject = setupObject.plotObject;
time = plotObject.time;

%formatting
a = plotObject.boundarySize*plotObject.deltaXY;
h = a*sqrt(3)/2;
lineWeight = plotObject.lineWeight;

%% getting nodal data
switch lower(time)
    case 'ref'
        nodes = obj.masterObject.qR;
    case 'n'
        nodes = obj.masterObject.qN;
    case 'n1'
        nodes = obj.masterObject.qN1;
    otherwise
        error('unknown timepoint')
end 

if isa(obj.masterObject, 'plateClass')
    nodes = obj.masterObject.meshObject.nodes;
end

%adujust dimension of arrays if 1D: everything in y-direction is zero
if obj.masterObject.dimension == 1 
    nodes = [nodes, zeros(size(nodes,1),1)];
end

%% plotting
nodesBC = obj.nodeList;
totBC = size(nodesBC,1);

x = zeros(totBC,4);
y = zeros(totBC,4);

for i = 1:numel(obj.nodalDof)
    switch obj.nodalDof(i)
        case 1
            xTri = [0, -h, -h, 0];
            yTri = [0, a/2, -a/2, 0];
        case 2
            xTri = [0, -a/2, a/2, 0];
            yTri = [0, -h, -h, 0];
        otherwise
            error('direction not implemented')
    end
    x(:,:) = kron(ones(totBC,1),xTri)+kron(ones(1,4),nodes(nodesBC,1));
    y(:,:) = kron(ones(totBC,1),yTri)+kron(ones(1,4),nodes(nodesBC,2));
    plot(x',y','r','linewidth',2*lineWeight)
end
end