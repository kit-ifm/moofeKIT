function ed = extract(edof,q)
%% Extrahieren der Elementdatentabelle aus der Elementfreiheitsgradzuordnungsdatentabelle und dem Geometrievektor q
%
% ed = extract(edof,q)
% 
% ed:   Elementdatentabelle
% edof: Elementfreiheitsgradzuordnungsdatentabelle
% q:    Geometrievektor

ed = q(edof')';

end