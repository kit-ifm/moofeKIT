function [nodes,edofTet,edofneumannTet] = linearTetrahedralCooksMembrane(anzElemX, anzElemY, anzElemZ, laenge, hoehelinks, hoeherechts, dicke)
%% Datei netz.m - Netzgenerator fuer triilineare isoparametrische 8-Knotenelemente
%   edof - Zuweisung der Knoten zum jeweiligen Element 
%   nodes - Lage der Knoten

%Verwendung 8 Knoten-Elementen
dim = 3;
anzNodes2D = (anzElemX+1)*(anzElemY+1);%nur eine Schicht (2D)
nodes = zeros(anzNodes2D*(anzElemZ+1),dim);
edof = zeros(anzElemX*anzElemY*anzElemZ,2^dim);
dx = laenge/anzElemX;
dz = dicke/anzElemZ;

%knotengeometrie ohne Riss
for i=1:anzElemX+1
   x=(i-1)*dx;
   hx=hoehelinks/laenge*x;
   dy=(hoeherechts/laenge*x+hoehelinks-hx)/anzElemY;
   
   nodes((i-1)*(anzElemY+1)+1:i*(anzElemY+1),1)=x;%
   nodes((i-1)*(anzElemY+1)+1:i*(anzElemY+1),2)=(0:anzElemY)*dy+hx;%y
end

for i=1:anzElemZ
    %x-y-Koordinate kopieren
    nodes(anzNodes2D*i+1:anzNodes2D*(i+1),1:2)=nodes(1:anzNodes2D,1:2);
    %z-Koordinate
    nodes(anzNodes2D*i+1:anzNodes2D*(i+1),3)=i*dz;
end

%edof 
for i=1:anzElemX
   edof((i-1)*anzElemY+1:i*anzElemY,1:4)= ...
       [(1:anzElemY)'+(i-1)*(anzElemY+1),(1:anzElemY)'+i*(anzElemY+1), ...
        (1:anzElemY)'+i*(anzElemY+1)+1,(1:anzElemY)'+(i-1)*(anzElemY+1)+1];
end

edof(1:anzElemX*anzElemY,5:8) = edof(1:anzElemX*anzElemY,1:4)+anzNodes2D;
for i=2:anzElemZ
    edof((i-1)*anzElemX*anzElemY+1:i*anzElemX*anzElemY,:) = edof(1:anzElemX*anzElemY,:)+anzNodes2D*(i-1);
end

%edofneumann
edofneumann=zeros(anzElemY*anzElemZ,4);

for i=1:anzElemZ
    edofneumann((i-1)*anzElemY+1:i*anzElemY,:)= ...
        [(1:anzElemY)'+(anzElemY+1)*anzElemX+(i-1)*anzNodes2D, ...
         (1:anzElemY)'+(anzElemY+1)*anzElemX+(i-1)*anzNodes2D+1, ...
         (1:anzElemY)'+(anzElemY+1)*anzElemX+(i)*anzNodes2D+1, ...
         (1:anzElemY)'+(anzElemY+1)*anzElemX+(i)*anzNodes2D];  
end

%--------------------------------------------------------------
% Einteilung eines Hexaeder-Elements in je 5 Tetraeder-Elemente
%--------------------------------------------------------------

newEl=size(edof,1)*5;
edofTet=zeros(newEl,4);
edofneumannTet=zeros(2*size(edofneumann,1),3);

for p=1:(newEl/5)
    edofTet(5*(p-1)+1,:)=edof(p,[8 5 1 6]);
    edofTet(5*(p-1)+2,:)=edof(p,[8 7 6 3]);
    edofTet(5*(p-1)+3,:)=edof(p,[3 4 1 8]);
    edofTet(5*(p-1)+4,:)=edof(p,[3 2 6 1]);
    edofTet(5*(p-1)+5,:)=edof(p,[3 6 8 1]);
end

for pp=1:size(edofneumann,1)
    edofneumannTet(2*(pp-1)+1,:)=edofneumann(pp,[1 2 4]);
    edofneumannTet(2*(pp-1)+2,:)=edofneumann(pp,[3 4 2]);
end

end
