function [nodes, edof] =  trilinearQuarterCylinder (r,h,nelPhi,nelR,nelZ)
    % Description
    % mesh for a 3D cylindrical body with inner Neumann domains
    % Philipp Kinon, 15.02.2018
    % based on circle-routines (meshCircle with corresponding functions) 
    %       by Robin Pfefferkorn


    % INPUT: 
    % r: radius (x-y-plane) h: height (z-direction)
    % nelPhi: number of elements in circumferential direction (one quadrant!)
    % nelR: number of elements in radial direction
    % nelZ: number of elements in z-direction

    % OUTPUT:
    % nodes: position of nodes 
    % edof: elements (trilinear; circle-routine also biquadratic possible but meshCyl only trilinear)

    %% Meshing one 2D circle
    [nodes2D,edof2D] = meshQuarterCircle(r,nelPhi,nelR,1,false);
    nel2D = size(edof2D,1);
    nno2D = size(nodes2D,1);

    %% Extruding the circle to a 3D cylinder
    %nodes
    nodesZDir = (0:h/nelZ:h)';
    nodes = [kron(ones(nelZ+1,1),nodes2D), kron(nodesZDir,ones(nno2D,1))];

    %edof
    edof = zeros(nel2D*nelZ,8);
    edof(1:nel2D,:) = [edof2D, edof2D+nno2D];
    edof(:,:) = kron(ones(nelZ,1),edof(1:nel2D,:))+kron((0:nelZ-1)',nno2D*ones(nel2D,8));

end
%###########################################################################################################################
%% SUBFUNCTIONS
%###########################################################################################################################
%% Function for meshing one quadrant of the 2D circle
function [nodes,edof] = meshQuarterCircle(r,nelPhi,nelR,order,useMiddlePart)
    %check input
    if mod(nelPhi,2)~=0
        error('nelPhi must be divedable by 2')
    end
    if useMiddlePart && mod(nelR,3)~=0 
        error('nelR must be divideable by 3')
    elseif ~useMiddlePart && mod(nelR,2)~=0
        error('nelR must be divideable by 2')
    end
    
    %geometry
    if useMiddlePart
        r1 = 2/3*r*1.0; %distance to inner border of outer circle
        r2 = 1/3*r*1.1; %distance to change to rectangle form
        r3 = r2*1.2; %distance to corner at 45deg of rectangle form
    else
        r = 2*r;
        r1 = 1/2*r; %outer border
        r2 = 1/2*r1;
        r3 = r2*1.2;
    end

    % number elements and nodes (total)
    % achtung: durch tie am schluss verringert sich die knotenzahl
    if ~useMiddlePart
        nelR=3/2*nelR;
    end
    nelRi = nelR/3; %elemete pro ring
    nel = nelPhi*2*nelRi+(nelPhi/2)^2;
    nno = (nelPhi*order+1)*(2*nelRi*order+1)+(nelPhi/2*order+1)^2;
    
    nodes = zeros(nno,2);
    edof = zeros(nel,(order+1)^2);

    phiMax=pi/2;
    
    % Aeusserer Ring (Kreisring) 
    %number of elements and nodes 
    nelX=nelPhi*(phiMax/(pi/2));
    nelY=nelRi;

    nnoX=nelX*order+1;
    nnoY=nelY*order+1;

    % Nodes
    dr=(r-r1)/nelRi/2^(order-1);
    dphi=(pi/2)/nelPhi/2^(order-1);
    phi=0:dphi:phiMax;
    ri=r:-dr:r1;

    x=cos(kron(ones(1,nnoY),phi)).*kron(ri,ones(1,nnoX));
    y=sin(kron(ones(1,nnoY),phi)).*kron(ri,ones(1,nnoX));

    nodes(1:nnoX*nnoY,1) = x;
    nodes(1:nnoX*nnoY,2) = y;
        
    % mittlerer Ring
    %points on inner border (edges)
    E=zeros(2*(phiMax/(pi/2))+1,2);
    ri=[r2, kron(ones(1,phiMax/(pi/2)),[r3,r2])];
    E(:,1)=cos(0:pi/4:phiMax).*ri;
    E(:,2)=sin(0:pi/4:phiMax).*ri;

    %all points on outer border
    A=zeros(nnoX,2);
    A(:,:)=nodes(nnoX*nnoY-nnoX+1:nnoX*nnoY,:);

    %all points on inner border
    nnoXsec=(nnoX-1)/2*(pi/2)/phiMax;

    B=zeros(nnoX,2);
    B(1,:)=E(1,:);
    for i=1:size(E,1)-1
        B(2+nnoXsec*(i-1):nnoXsec*i+1,:)=createPoints(E(i,:),E(i+1,:),nnoXsec,false);
    end

    %nodes
    help=createPoints(A,B,nnoY-1,false);

    nodes(nnoX*nnoY+1:nnoX*(2*nnoY-1),1)=help(:,1);
    nodes(nnoX*nnoY+1:nnoX*(2*nnoY-1),2)=help(:,2);

    %edof
    edof(1:nelX*2*nelY,:) = meshEdofRectangle(1,nelX,2*nelY,order);

    % Mittelteil (rechteck)
    node1=nnoX*(2*nnoY-1)+1;
    elem1=nelX*2*nelY+1;

    nelX = nelPhi/2;
    nelY = nelPhi/2;
    nnoX = nelX*order+1;
    nnoY = nelY*order+1;

    %points on right border (taken from middle part)
    A=B(1:nnoX,:);

    %points on left border
    B=zeros(nnoX,2);

    B(:,1) = 0;
    B(:,2) = 0:r2/(nnoX-1):r2;

    %nodes
    help = createPoints(A,B,nnoY,true);
    nodes(node1:node1+nnoX*nnoY-1,:) = help;

    %edof
    edof(elem1:elem1+nelX*nelY-1,:) = meshEdofRectangle(node1,nelX,nelY,order);
    
    % delete surplus elements and nodes
    if ~useMiddlePart
        nelX=nelPhi*(phiMax/(pi/2));
        nelY=nelRi;

        nnoX=nelX*order+1;
        nnoY=nelY*order+1;
        
        nodes(1:nnoX*nnoY-nnoX,:) = [];
        edof(1:nelX*nelY,:) = [];
        edof = edof-(nnoX*nnoY-nnoX);
    end

    % Tie mesh parts (inner part not connected to outer part so far)
    [nodes, edof] = meshTie(nodes,edof,zeros(0,size(nodes,2)),zeros(0,size(edof,2)),r/1000);
end

%% Function for creating individual points
function [ points ] = createPoints(A,B,numPoints,includeA)
    %pointsOnLine calculates points on a line
    %   A...first point on line
    %   B...second point on line
    %   numPoints...number of points along line
    %   includA...if A is also included

    if size(A,1)~=size(B,1)
        error('A und B muessen die gleiche Dimension haben')
    end

    if includeA 
        lam = 0:1/(numPoints-1):1;
        points = zeros(numPoints*size(A,1),2);
        points(:,:) = kron(lam',B-A)+kron(ones(1,numPoints)',A);
    else
        lam = 1/numPoints:1/numPoints:1;
        points = zeros(numPoints*size(A,1),2);
        points(:,:) = kron(lam',B-A)+kron(ones(1,numPoints)',A);
    end
end

%% Function for meshing the rectangular parts in the quadrant
function [edof] = meshEdofRectangle(node1,nelX,nelY,order)
    % number elements and nodes (total)
    nel=nelX*nelY;
    nnoX=nelX*order+1;
    nnoY=nelY*order+1;

    % edof
    if order == 1 %bilinear shape functions    
        edof = zeros(nel,4);
        edof(:,1)=kron(ones(nelY,1),(1:(nnoX-1))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelX,1));
        edof(:,2)=kron(ones(nelY,1),(2:(nnoX-0))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelX,1));
        edof(:,3)=kron(ones(nelY,1),(2:(nnoX-0))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelX,1));
        edof(:,4)=kron(ones(nelY,1),(1:(nnoX-1))')+kron((0:nnoX:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelX,1));
    elseif order == 2 %biquadratic shape functions
        edof = zeros(nel,9);
        edof(:,1)=kron(ones(nelY,1),(1:2:(nnoX-2))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelX,1));
        edof(:,2)=kron(ones(nelY,1),(3:2:(nnoX-0))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelX,1));
        edof(:,3)=kron(ones(nelY,1),(3:2:(nnoX-0))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+2*nnoX,ones(nelX,1));
        edof(:,4)=kron(ones(nelY,1),(1:2:(nnoX-2))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+2*nnoX,ones(nelX,1));
        edof(:,5)=kron(ones(nelY,1),(2:2:(nnoX-1))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+0*nnoX,ones(nelX,1));
        edof(:,6)=kron(ones(nelY,1),(3:2:(nnoX-0))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelX,1));
        edof(:,7)=kron(ones(nelY,1),(2:2:(nnoX-1))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+2*nnoX,ones(nelX,1));
        edof(:,8)=kron(ones(nelY,1),(1:2:(nnoX-2))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelX,1));
        edof(:,9)=kron(ones(nelY,1),(2:2:(nnoX-1))')+kron((0:nnoX*2:(nnoX-1)*(nnoY-1))'+1*nnoX,ones(nelX,1));
    else
        error('unable to mesh: elementtype currently not implemented') 
    end

    edof = edof+node1-1;
end
