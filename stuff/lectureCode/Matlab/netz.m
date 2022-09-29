function [edof,q] = netz(elementx, elementy, laenge, hoehe)
%% Datei netz.m - Netzgenerator fuer bilineare isoparametrische 4-Knotenelemente
%   edof - Zuweisung der Knoten zum jeweiligen Element 
%   q    - Lage der Knoten
%
%
%     (5,6) o--------------o- - -
%           |              |            
%           |              |   
%           |     (2)      |      
%           |              |              |              
%           |              | (9,10)       |              
%     (3,4) o--------------o--------------o- - -
%           |              |              |     
%           |              |              |     
%           |     (1)      |      (3)     |    
%           |              |              |       
%           |              |              |             
%           o--------------o--------------o- - -
%     (1,2)              (7,8)



nel      = elementx*elementy;    % Anzahl Gesamtelemente 
k        = (elementy+1)*2;       % Anzahl Knotenfreiheitsgrade in y-Richtung
hx       = laenge/elementx;      % Elementlaenge in x-Richtung

% Allozieren
edof = zeros(nel,8);
ex   = zeros(nel,4);
ey   = zeros(nel,4);

% Zuweisung Knotenfreiheitsgrade <-> Element
for ii = 1:elementx
    edof((1:elementy)+(ii-1)*elementy,[1 2 7 8]) = (((1:elementy)'-1)*2 + (ii-1)*k) * ones(1,4)   + ones(elementy,1)*(1:4);
    edof((1:elementy)+(ii-1)*elementy,[3 4 5 6]) = (((1:elementy)'-1)*2 + (ii-1)*k+k) * ones(1,4) + ones(elementy,1)*(1:4);
end

% Lage der Knoten in x - Richtung
for jj = 1:elementx
    ex((1:elementy)+(jj-1)*elementy,1:4) = [hx*((ones(elementy,1))*jj-1) hx*((ones(elementy,1))*jj) hx*((ones(elementy,1))*jj) hx*((ones(elementy,1))*jj-1)];
end

% Lage der Knoten in y - Richtung
laenge_element = laenge/elementx;
Y1             = laenge_element*(hoehe/laenge);
Y2             = laenge_element*(16/laenge);

for l = 1:elementx
    y1 = Y1*(l-1);
    y2 = Y1*l;
    y3 = Y2*(l-1) + hoehe;
    y4 = Y2*l + hoehe;
    
    ya = (y3-y1)/elementy;
    yb = (y4-y2)/elementy;
    
    for m = 1:elementy
        ey(m+(l-1)*elementy,1) = y1 + ya*(m-1);
        ey(m+(l-1)*elementy,2) = y2 + yb*(m-1);
        ey(m+(l-1)*elementy,3) = y2 + yb*m;
        ey(m+(l-1)*elementy,4) = y1 + ya*m;
    end
end

ed = [ex(:,1) ey(:,1) ex(:,2) ey(:,2) ex(:,3) ey(:,3) ex(:,4) ey(:,4)];


% Geometrievektor erstellen
[nie,n]=size(edof);
t = edof(:,1:n);
q = zeros(2*(elementx+1)*(elementy+1),1);
for i = 1:nie
    q(t(i,:),1) = ed(i,1:(n))';
end

