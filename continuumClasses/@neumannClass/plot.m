function plot(obj,setupObject)
%plots loads of neumann
plotObject = setupObject.plotObject;
time = plotObject.time;

updateNodalForce(obj)

%formatting
lineWeight = plotObject.lineWeight;
deltaXY = plotObject.deltaXY;
lengthPoint = 0.2;
phiPoint = pi/6;

%% getting nodal data
nodalForce = obj.nodalForce;
maxF = max(max(abs(nodalForce)));
switch time
    case 'R'
        nodes = obj.masterObject.qR;
    case 'N'
        nodes = obj.masterObject.qN;
    case 'N1'
        nodes = obj.masterObject.qN1;
    otherwise
        error('unknown timepoint')
end 

%adujust dimension of arrays if 1D: everything in y-direction is zero
if obj.masterObject.dimension==1 
    nodalForce = [nodalForce, zeros(size(nodalForce,1),1)];
    nodes = [nodes, zeros(size(nodes,1),1)];
end

numLoad = size(nodalForce,1);

%% plotting of mesh
x = zeros(numLoad,5);
y = zeros(numLoad,5);

for e=1:numLoad    
    Fx = nodalForce(e,1);
    Fy = nodalForce(e,2);
    F = sqrt(Fx^2+Fy^2);
    
    if F>eps
        points = zeros(2,5);
        L = F/maxF*0.08*deltaXY;
        
        points(1,1)= L;
        points(1,2)= 0;
        points(1,3)= L*lengthPoint;
        points(1,4)= 0;
        points(1,5)= L*lengthPoint;

        points(2,1)= 0;
        points(2,2)= 0;
        points(2,3)= L*lengthPoint*tan(phiPoint);
        points(2,4)= 0;
        points(2,5)= -L*lengthPoint*tan(phiPoint);
        
        %rotation and translation
        sinP = Fy/F;
        cosP = -Fx/F;
        A = [cosP, sinP; -sinP, cosP];
        points = A*points;
        
        x(e,:)=points(1,:)+nodes(obj.masterNodes(e),1);
        y(e,:)=points(2,:)+nodes(obj.masterNodes(e),2);
    end
end

plot(x',y','b','linewidth',lineWeight)
end