function [ out ] = findNodes(nodes, line, tol)
%gives indices of nodes that within tol to line
%   input: nodes (all nodes from which subset should be selected)
%   input: line (2D-array with points on which nodes are searched)
%   output: indices of nodes from nodes on line 

%% checks
if size(line,1)<2 
    error('line needs at least two points')
end
if size(line,2)~=2
    error('line must have have two columns (x and y)')
end
if size(nodes,2)~=2
    error('nodes must have have two columns (x and y)')
end
if tol<0
    error('tol must be greater than zero')
end

%% find nodes
%general values
nno = size(nodes,1);
onLine = logical(ones(nno,1)*false);

% loop over all segments of line
for i=1:size(line,1)-1
    x1 = line(i,:)';
    x2 = line(i+1,:)';
    
    %vector of linesegment and perpendicular vector
    r = x2-x1;
    p = [-r(2);r(1)];
    
    A = [r,-p];
    
    %calculate intersections
    lam = A\(nodes'-kron(ones(1,nno),x1));
    lam(2,:) = [];
    
    %determine if nodes within tolerance (three cases)
    %case one: intersection before x1
    act = (lam<=0);
    no1 = nodes(act,:);
    no2 = kron(ones(size(no1,1),1),x1');
    onLine(act) = onLine(act) | (sum((no1-no2).^2,2)<=tol);
    %case two: intersection after x2
    act = (lam>=1);
    no1 = nodes(act,:);
    no2 = kron(ones(size(no1,1),1),x2');
    onLine(act) = onLine(act) | (sum((no1-no2).^2,2)<=tol);
    %case three: inbetween
    act = (lam>0 & lam<1);
    no1 = nodes(act,:);
    no2 = kron(ones(size(no1,1),1),x1')+kron(lam(act)',r');
    onLine(act) = onLine(act) | (sum((no1-no2).^2,2)<=tol);
end

ind = (1:nno)';
out = ind(onLine);

end

