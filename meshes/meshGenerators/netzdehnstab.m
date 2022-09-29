function [nodes,edof,edofneumann] = netzdehnstab(anzElemX, laenge)
% %% Datei netz.m - Netzgenerator fuer 1D-Dehnstab
% %   edof - Zuweisung der Knoten zum jeweiligen Element 
% %   nodes - Lage der Knoten
% 
% %Verwendung 2 Knoten-Elementen
% % dim = 1;
% % anzNodes2D = (anzElemX+1)*(anzElemY+1);%nur eine Schicht (2D)
% % nodes = zeros(anzNodes2D*(anzElemZ+1),dim);
% edof = zeros(anzElemX,2);
% dx = laenge/anzElemX;
% 
% %knotengeometrie ohne Riss
% nodes = zeros((anzElemX+1),1);
% nodes(1:size(nodes,1))=(0:anzElemX)*dx;
% 
% %edof 
% for i=1:anzElemX
%    edof(i,:)= [i i+1];
% end
% 
% %edofneumann
% edofneumann=anzElemX+1;

warning('netzdehnstab() is deprecated and will be removed soon. Please consider using meshOneDimensional() instead.');

[nodes, edof, edofneumann] = meshOneDimensional(laenge, anzElemX, 1);
nodes = nodes + laenge/2;

end



