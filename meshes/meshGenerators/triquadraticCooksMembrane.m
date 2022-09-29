function [nodes,edof,edofneumann] = triquadraticCooksMembrane(anzElemX, anzElemY, anzElemZ, laenge, hoehelinks, hoeherechts, dicke)

%%Startscript%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functionname:     cookQ23D.m
% Creation date:    08.12.2016
% Creator:          P. Kinon 
% Discription:      3D-mesh-generation for 27-node brick elements in a
%                   hexahedron geometry (e.g. cooks membrane).
%                   Edof, nodes and neumann-edof are created by using
%                   geometry data.
% Based on:         2D-mesh-generation-function provided in "Grundlagen FE"

% Modifications:    --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nel = anzElemX*anzElemY*anzElemZ; % Anzahl Elemente
k0 = (anzElemX-1)*6+9;
k1 = (anzElemY-1)*6+9;  % Anzahl an FHG'en der x-y-Ebene auf Gerade x=const
k2 = k1*(2*anzElemZ+1); % Anzahl aller FHG'e auf der Ebene x=const
nFHG = k2*(2*anzElemX+1); % Anzahl aller FHG
nodesy=2*anzElemY+1;
nodesx=2*anzElemX+1;
nedof=zeros(nel,27*3);
nedof(1:nel,1) = (1:nel)'; % Durchnummerieren der Elemente in der ersten Spalte

%---------------------------------------------------------------
% Zuweisung Knoten <-> Element, edof
%---------------------------------------------------------------

% Hilfsmatrix zur Knotendefinition

nmbrs=zeros(nodesy,3*nodesx);
nmbrs=ones(nodesy,1)*(1:3)+(0:3:3*(nodesy-1))'*[1 1 1];
nmbrstemp=nmbrs;

for kk=1:(nodesx-1)
    nmbrs=[nmbrs nmbrstemp+3*kk*nodesy];
end

nmbrs_temp=nmbrs;
for ll=1:anzElemZ
    nmbrs=[nmbrs; nmbrs_temp+k0*(2*ll-1)*nodesy; nmbrs_temp+2*k0*ll*nodesy];
end

% Addressierung der FHG-->nedof

for i=1:anzElemZ
    for j=1:anzElemY
        for k=1:anzElemX
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[1 2 3 25 26 27 13 14 15]) = nmbrs( (2*j-1)+nodesy*(i-1)*2, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[37 38 39 67 68 69 52 53 54]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+1, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[4 5 6 28 29 30 16 17 18]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+2, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[46 47 48 76 77 78 61 62 63]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+nodesy, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[49 50 51 79 80 81 64 65 66]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+nodesy+1, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[40 41 42 70 71 72 55 56 57]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+nodesy+2, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[10 11 12 34 35 36 22 23 24]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+2*nodesy, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[43 44 45 73 74 75 58 59 60]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+2*nodesy+1, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
            nedof( (i-1)*anzElemY*anzElemX+(j-1)*anzElemX+k ,[7 8 9 31 32 33 19 20 21]) = nmbrs( (2*j-1)+nodesy*(i-1)*2+2*nodesy+2, [1 2 3 4 5 6 7 8 9]+(k-1)*6);
        end
    end
end

edof=nedof(:,(3:3:size(nedof,2)))/3;

%---------------------------------------------------------------
% Lage der Knoten, q-Vektor nodes 
%---------------------------------------------------------------

Nelx=anzElemX*2; % Anzahl Knotenzwischenraeume in x-Richtung
Nely=anzElemY*2; % Anzahl Knotenzwischenraeume in y-Richtung

Cx=[0 laenge laenge 0];  % Positionen x-Richtung der Eckknoten 2D
Cy=[0 hoehelinks hoeherechts hoehelinks]; % Positionen y-Richtung der Eckknoten 2D

% Erzeugen der Knotenkoordinaten

e=1;
for ii=1:Nelx+1
    for j=1:Nely+1
        xQ(e)=Cx(1)+(Cx(4)-Cx(1))*(j-1)/Nely+(Cx(2)-Cx(1)+(Cx(1)+Cx(3)-Cx(2)-Cx(4))*(j-1)/Nely)*(ii-1)/Nelx; % x-Koordinaten
        yQ(e)=Cy(1)+(Cy(4)-Cy(1))*(j-1)/Nely+(Cy(2)-Cy(1)+(Cy(1)+Cy(3)-Cy(2)-Cy(4))*(j-1)/Nely)*(ii-1)/Nelx; % y-Koordinaten
        e=e+1;
    end
end

xQQ=xQ;
yQQ=yQ;
zQ=zeros(1,size(xQ,2));
zQ0=zQ;

for jj=1:anzElemZ
xQ=[xQ xQQ xQQ];
yQ=[yQ yQQ yQQ];
zQ=[zQ (2*jj-1)*dicke/(2*anzElemZ)+zQ0 jj*dicke/anzElemZ+zQ0];
end

nodes=[xQ' yQ' zQ'];

%---------------------------------------------------------------
% Neumann-Edof
%---------------------------------------------------------------

edofneumann=zeros(anzElemY*anzElemZ,9);

g=1;
for i=1:anzElemZ
    for j=1:anzElemY
       edofneumann(g,[1 8 4])= (nodesx-1)*nodesy+[1 2 3]+(j-1)*2+(i-1)*2*nodesy*nodesx;
       edofneumann(g,[5 9 7])= (nodesx-1)*nodesy+[1 2 3]+(j-1)*2+(i-1)*2*nodesy*nodesx+nodesx*nodesy;
       edofneumann(g,[2 6 3])= (nodesx-1)*nodesy+[1 2 3]+(j-1)*2+(i-1)*2*nodesy*nodesx+2*nodesx*nodesy;
       g=g+1;
    end
end

end